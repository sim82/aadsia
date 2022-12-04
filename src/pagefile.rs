use rand::Rng;
use std::{
    fs::{File, OpenOptions},
    io::{Read, Seek, SeekFrom, Write},
    os::unix::prelude::FileExt,
    path::Path,
};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum Error {
    #[error("underlying io error")]
    IoError(#[from] std::io::Error),
}

pub type Result<T> = std::result::Result<T, Error>;

pub trait PageFile {
    fn alloc_page(&mut self) -> u64;
    fn free_page(&mut self, id: u64);
    fn read_page(&mut self, id: u64) -> Result<Vec<u8>>;
    fn write_page(&mut self, id: u64, data: &[u8]);
}

const ENT_TERM: u64 = u64::MAX;

pub struct PageFileStd {
    file: File,
    page_size: usize,
    last_page: u64,
}

impl PageFileStd {
    pub fn open<P: AsRef<Path>>(filename: P, page_size: usize) -> Result<Self> {
        let file = OpenOptions::new().read(true).write(true).open(filename)?;

        let metadata = file.metadata()?;
        let file_size = metadata.len();
        assert!(file_size % page_size as u64 == 0);
        assert!(file_size >= page_size as u64);
        let last_page = file_size / page_size as u64 - 1;
        println!("last page: {}", last_page);

        Ok(Self {
            file,
            page_size,
            last_page,
        })
    }
    pub fn create<P: AsRef<Path>>(filename: P, page_size: usize) -> Result<Self> {
        let file = OpenOptions::new()
            .read(true)
            .write(true)
            .truncate(true)
            .create(true)
            .open(filename)?;

        file.set_len(page_size as u64);

        let mut almost_done = Self {
            file,
            page_size,
            last_page: 0,
        };
        almost_done.write_freelist(0, ENT_TERM);
        Ok(almost_done)
    }

    fn read_freelist(&mut self, id: u64) -> u64 {
        let mut buf = [0u8; 8];

        self.file
            .read_exact_at(&mut buf, self.page_offset(id) as u64)
            .unwrap();

        // self.file.seek(SeekFrom::Start(self.page_offset(id) as u64));
        // self.file.read_exact(&mut buf);

        u64::from_le_bytes(buf)
    }

    fn write_freelist(&mut self, id: u64, next_id: u64) {
        let mut buf = [0u8; 8];
        self.file
            .write_all_at(&u64::to_le_bytes(next_id), self.page_offset(id) as u64);

        // self.file.seek(SeekFrom::Start(self.page_offset(id) as u64));
        // self.file.write_all(&u64::to_le_bytes(next_id));
    }

    fn page_offset(&self, id: u64) -> usize {
        self.page_size * id as usize
    }

    fn dump_freelist(&mut self) {
        let mut next = 0;

        while next != ENT_TERM {
            next = self.read_freelist(next);
            println!("{:x}", next);
        }
    }
}

impl PageFile for PageFileStd {
    fn alloc_page(&mut self) -> u64 {
        let next_free = self.read_freelist(0);

        if next_free != ENT_TERM {
            let next_next_free = self.read_freelist(next_free);
            // let page = self.freelist_page(next_next_free);
            // self.write_page(0, &page);
            self.write_freelist(0, next_next_free);

            return next_free;
        } else {
            // self.file.
            self.last_page += 1;
            self.file
                .set_len((self.last_page + 1) * self.page_size as u64);

            return self.last_page;
        }
    }

    fn free_page(&mut self, id: u64) {
        let head_next = self.read_freelist(0);
        // let page = self.freelist_page(head_next);
        // self.write_page(id, &page);
        // let page = self.freelist_page(id);
        // self.write_page(0, &page);
        self.write_freelist(id, head_next);
        self.write_freelist(0, id);
    }

    fn read_page(&mut self, id: u64) -> Result<Vec<u8>> {
        todo!()
    }

    fn write_page(&mut self, id: u64, data: &[u8]) {
        self.file
            .write_at(data, self.page_offset(id) as u64)
            .unwrap();
    }
}

#[test]
fn test_create() {
    let mut pf = PageFileStd::create("/tmp/test.bin", 16).unwrap();
    let page1 = pf.alloc_page();
    let page2 = pf.alloc_page();

    println!("page1: {}", page1);
    println!("page2: {}", page2);

    pf.free_page(page1);
    pf.free_page(page2);

    let page1 = pf.alloc_page();
    let page2 = pf.alloc_page();
    let page3 = pf.alloc_page();

    println!("page1: {}", page1);
    println!("page2: {}", page2);
    println!("page3: {}", page3);

    pf.free_page(page1);
    pf.free_page(page2);
    pf.free_page(page3);
}

#[test]
fn test_open() {
    let mut pf = PageFileStd::open("/tmp/test.bin", 16).unwrap();
    pf.dump_freelist();
    let page1 = pf.alloc_page();
    let page2 = pf.alloc_page();

    println!("page1: {}", page1);
    println!("page2: {}", page2);
    pf.dump_freelist();
}

#[test]
fn test_random() {
    let mut pf = PageFileStd::create("/tmp/test.bin", 16).unwrap();

    let mut pages = Vec::new();

    let mut rng = rand::thread_rng();

    for _ in 0..100000 {
        if !pages.is_empty() && rng.gen_bool(0.5) {
            let delete_index = rng.gen_range(0..pages.len());

            let last_index = pages.len() - 1;
            if delete_index < last_index {
                pages.swap(delete_index, last_index);
            }
            let to_delete = pages.pop().unwrap();
            println!("free: {}", to_delete);
            pf.free_page(to_delete as u64);
        } else {
            let page = pf.alloc_page();
            println!("alloc: {}", page);
            pages.push(page);
        }
    }
}
