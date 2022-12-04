use aadsia::pagefile::{PageFile, PageFileStd};
use rand::Rng;

fn main() {
    let mut pf = PageFileStd::create("/tmp/test.bin", 16).unwrap();

    let mut pages = Vec::new();

    let mut rng = rand::thread_rng();

    for _ in 0..1000000 {
        if !pages.is_empty() && rng.gen_bool(0.3) {
            let delete_index = rng.gen_range(0..pages.len());

            let last_index = pages.len() - 1;
            if delete_index < last_index {
                pages.swap(delete_index, last_index);
            }
            let to_delete = pages.pop().unwrap();
            // println!("free: {}", to_delete);
            pf.free_page(to_delete as u64);
        } else {
            let page = pf.alloc_page();
            // println!("alloc: {}", page);
            pages.push(page);
        }
    }
}
