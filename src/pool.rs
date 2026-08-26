//! A tiny ordered worker pool.
//!
//! Jobs are submitted from a single producer, executed on `n` worker threads and
//! handed back to the consumer **in submission order**. Both the BGZF reader and
//! the BGZF writer are built on top of it: block (de)compression is embarrassingly
//! parallel but the byte stream it belongs to is strictly ordered.

use std::collections::HashMap;
use std::sync::mpsc::{sync_channel, Receiver, SyncSender};
use std::sync::{Arc, Condvar, Mutex};
use std::thread::JoinHandle;

struct Slots<O> {
    done: HashMap<u64, O>,
    /// Set once every worker has exited, so a consumer waiting on a sequence
    /// number that will never arrive can bail out instead of blocking forever.
    workers_live: usize,
}

pub struct OrderedPool<O: Send + 'static> {
    tx: Option<SyncSender<(u64, Box<dyn FnOnce() -> O + Send>)>>,
    slots: Arc<(Mutex<Slots<O>>, Condvar)>,
    handles: Vec<JoinHandle<()>>,
    next_in: u64,
    next_out: u64,
    /// Submitted but not yet consumed.
    outstanding: usize,
    capacity: usize,
}

impl<O: Send + 'static> OrderedPool<O> {
    /// `threads` workers, allowing up to `depth` jobs in flight.
    pub fn new(threads: usize, depth: usize) -> Self {
        let threads = threads.max(1);
        let capacity = depth.max(threads * 2);
        let (tx, rx) = sync_channel::<(u64, Box<dyn FnOnce() -> O + Send>)>(capacity);
        let rx = Arc::new(Mutex::new(rx));
        let slots = Arc::new((
            Mutex::new(Slots { done: HashMap::new(), workers_live: threads }),
            Condvar::new(),
        ));

        let mut handles = Vec::with_capacity(threads);
        for _ in 0..threads {
            let rx: Arc<Mutex<Receiver<_>>> = Arc::clone(&rx);
            let slots = Arc::clone(&slots);
            handles.push(std::thread::spawn(move || {
                loop {
                    // Hold the receiver lock only long enough to take one job.
                    let job = { rx.lock().unwrap().recv() };
                    match job {
                        Ok((seq, f)) => {
                            let out = f();
                            let (m, cv) = &*slots;
                            m.lock().unwrap().done.insert(seq, out);
                            cv.notify_all();
                        }
                        Err(_) => break,
                    }
                }
                let (m, cv) = &*slots;
                m.lock().unwrap().workers_live -= 1;
                cv.notify_all();
            }));
        }

        OrderedPool {
            tx: Some(tx),
            slots,
            handles,
            next_in: 0,
            next_out: 0,
            outstanding: 0,
            capacity,
        }
    }

    pub fn capacity(&self) -> usize {
        self.capacity
    }

    pub fn outstanding(&self) -> usize {
        self.outstanding
    }

    /// Queue a job. Blocks once `capacity` jobs are in flight, which is what
    /// keeps memory bounded on multi-gigabyte inputs.
    pub fn submit<F: FnOnce() -> O + Send + 'static>(&mut self, f: F) {
        let seq = self.next_in;
        self.next_in += 1;
        self.outstanding += 1;
        self.tx
            .as_ref()
            .expect("pool closed")
            .send((seq, Box::new(f)))
            .expect("worker threads died");
    }

    /// Take the next result in submission order, or `None` when everything
    /// submitted so far has been consumed.
    pub fn next(&mut self) -> Option<O> {
        if self.outstanding == 0 {
            return None;
        }
        let (m, cv) = &*self.slots;
        let mut g = m.lock().unwrap();
        loop {
            if let Some(v) = g.done.remove(&self.next_out) {
                self.next_out += 1;
                self.outstanding -= 1;
                return Some(v);
            }
            if g.workers_live == 0 {
                return None;
            }
            g = cv.wait(g).unwrap();
        }
    }
}

impl<O: Send + 'static> Drop for OrderedPool<O> {
    fn drop(&mut self) {
        self.tx.take();
        for h in self.handles.drain(..) {
            let _ = h.join();
        }
    }
}
