use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use rust_htslib::bam::record::{Record, Aux};

use rust_htslib::bam::Read;
use itertools::Itertools;

fn main() {
    let args: Vec<String> = env::args().collect();
    let bamfile = &args[1];
    let outfolder = &args[2];
    let barcode_tag = &args[3];
    let bcfile = &args[4];

    let base = Path::new(bamfile).file_stem().unwrap().to_str().unwrap();
    let basename = base.to_string();

    fn get_barcode(read: &Record, barcode_tag: &str) -> String {
        match read.aux(barcode_tag.as_bytes()) {
            Ok(tag) => if let Aux::String(v) = tag {v.to_string()} else {"NA".to_string()},
            Err(_e) => "NA".to_string()
        }
    }

    fn write_passing_reads(bc_dict: &fnv::FnvHashMap<String, usize>, bamfile: &str, files: &mut Vec<rust_htslib::bam::Writer>, barcode_tag: &str) {
        let mut bam = rust_htslib::bam::Reader::from_path(&bamfile).unwrap();
        let Itr = bam.records();

        for record in Itr {
            let read = record.unwrap();
            let read_barcode = get_barcode(&read, barcode_tag);

            // If read barcode is in whitelist, then write it out
            if bc_dict.contains_key(&read_barcode) {
                let idx = bc_dict[&read_barcode];
                let file = &mut files[idx];
                file.write(&read).unwrap();
            }
        }
    }

    // Read in the barcodes
    let barcode_file_handle = File::open(bcfile).unwrap();
    let barcode_file_reader = BufReader::new(barcode_file_handle);
    let bc: Vec<String> = barcode_file_reader.lines().map(|x| x.unwrap().trim().to_owned()).collect();

    // Open up a bunch of files and write out reads for valid barcodes
    let bambcfiles: Vec<String> = bc.iter().map(|bc1| format!("{}/{}.bam", outfolder, bc1)).collect();
    let bc_dict: fnv::FnvHashMap<String, usize> = bc.iter().enumerate().map(|(i, bc1)| (bc1.to_owned(), i)).collect();
    let _temp = rust_htslib::bam::Reader::from_path(&bamfile).unwrap();
    let mut files = bambcfiles.iter().map(|filename| rust_htslib::bam::Writer::from_path(filename, &rust_htslib::bam::header::Header::from_template(&_temp.header().clone()), rust_htslib::bam::Format::Bam).unwrap()).collect_vec();
    drop(_temp);

    // Final loop to write out passing reads
    write_passing_reads(&bc_dict, &bamfile, &mut files, &barcode_tag);
}