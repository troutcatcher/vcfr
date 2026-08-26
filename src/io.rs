//! Input/output plumbing: format sniffing, threaded codecs, stdin/stdout.

use std::fs::File;
use std::io::{self, BufReader, BufWriter, Cursor, Read, Write};

use crate::bgzf::reader::{BgzfReader, SerialBgzfReader};
use crate::bgzf::{is_bgzf_header, is_gzip_header, BgzfWriter, RawBlockReader};
use crate::lines::{BufSource, LineReader, PlainSource};

impl BufSource for Box<dyn BufSource + Send> {
    fn next_chunk(&mut self) -> io::Result<Option<Vec<u8>>> {
        (**self).next_chunk()
    }
    fn recycle(&mut self, buf: Vec<u8>) {
        (**self).recycle(buf)
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Compression {
    None,
    Bgzf,
    Gzip,
}

pub type Reader = LineReader<Box<dyn BufSource + Send>>;

/// Reads enough of `r` to identify the compression, then splices those bytes
/// back onto the front of the stream.
fn sniff<R: Read + Send + 'static>(mut r: R) -> io::Result<(Compression, Box<dyn Read + Send>)> {
    let mut head = vec![0u8; 64];
    let mut n = 0;
    while n < head.len() {
        match r.read(&mut head[n..])? {
            0 => break,
            k => n += k,
        }
    }
    head.truncate(n);
    let comp = if is_bgzf_header(&head) {
        Compression::Bgzf
    } else if is_gzip_header(&head) {
        Compression::Gzip
    } else {
        Compression::None
    };
    Ok((comp, Box::new(Cursor::new(head).chain(r))))
}

fn raw_source(path: &str) -> io::Result<Box<dyn Read + Send>> {
    if path == "-" {
        Ok(Box::new(io::stdin()))
    } else {
        Ok(Box::new(BufReader::with_capacity(
            1 << 20,
            File::open(path).map_err(|e| io::Error::new(e.kind(), format!("{path}: {e}")))?,
        )))
    }
}

/// Open a VCF for reading, transparently handling plain, gzip and BGZF input.
///
/// `workers` BGZF inflation threads are spawned; `0` inflates inline on the
/// calling thread. Even a single worker is worth having when the caller has
/// other work to do, since it decouples inflation from consumption.
pub fn open_reader(path: &str, workers: usize) -> io::Result<Reader> {
    let (comp, r) = sniff(raw_source(path)?)?;
    let src: Box<dyn BufSource + Send> = match comp {
        Compression::Bgzf if workers >= 1 => Box::new(BgzfReader::new(r, workers)),
        Compression::Bgzf => Box::new(SerialBgzfReader::new(r)),
        Compression::Gzip => Box::new(PlainSource::new(flate2::read::MultiGzDecoder::new(r))),
        Compression::None => Box::new(PlainSource::new(r)),
    };
    Ok(LineReader::new(src))
}

pub fn detect_compression(path: &str) -> io::Result<Compression> {
    let (c, _) = sniff(raw_source(path)?)?;
    Ok(c)
}

/// Open a BGZF file for verbatim block access (`concat --naive`).
pub fn open_raw_blocks(path: &str) -> io::Result<RawBlockReader<Box<dyn Read + Send>>> {
    let (comp, r) = sniff(raw_source(path)?)?;
    if comp != Compression::Bgzf {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("{path}: not BGZF-compressed"),
        ));
    }
    Ok(RawBlockReader::new(r))
}

/// Output sink: either plain bytes or a threaded BGZF stream.
pub enum Writer {
    Plain(BufWriter<Box<dyn Write + Send>>),
    Bgzf(BgzfWriter<BufWriter<Box<dyn Write + Send>>>),
}

impl Writer {
    pub fn create(path: Option<&str>, bgzf: bool, threads: usize, level: u32) -> io::Result<Writer> {
        let sink: Box<dyn Write + Send> = match path {
            None | Some("-") => Box::new(io::stdout()),
            Some(p) => Box::new(
                File::create(p).map_err(|e| io::Error::new(e.kind(), format!("{p}: {e}")))?,
            ),
        };
        let buffered = BufWriter::with_capacity(1 << 20, sink);
        Ok(if bgzf {
            Writer::Bgzf(BgzfWriter::new(buffered, threads, level))
        } else {
            Writer::Plain(buffered)
        })
    }

    /// Append a raw BGZF block. Only valid on a BGZF sink.
    pub fn write_raw_block(&mut self, block: Vec<u8>) -> io::Result<()> {
        match self {
            Writer::Bgzf(w) => w.write_raw_block(block),
            Writer::Plain(_) => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "raw block passthrough requires BGZF output",
            )),
        }
    }

    pub fn finish(&mut self) -> io::Result<()> {
        match self {
            Writer::Plain(w) => w.flush(),
            Writer::Bgzf(w) => w.finish(),
        }
    }
}

impl Write for Writer {
    #[inline]
    fn write(&mut self, data: &[u8]) -> io::Result<usize> {
        match self {
            Writer::Plain(w) => w.write(data),
            Writer::Bgzf(w) => w.write(data),
        }
    }
    #[inline]
    fn write_all(&mut self, data: &[u8]) -> io::Result<()> {
        match self {
            Writer::Plain(w) => w.write_all(data),
            Writer::Bgzf(w) => w.write_all(data),
        }
    }
    fn flush(&mut self) -> io::Result<()> {
        match self {
            Writer::Plain(w) => w.flush(),
            Writer::Bgzf(w) => w.flush(),
        }
    }
}

/// Whether a path implies BGZF output, used when `--output-type` is left unset.
pub fn path_implies_bgzf(path: Option<&str>) -> bool {
    matches!(path, Some(p) if p.ends_with(".gz") || p.ends_with(".bgz"))
}
