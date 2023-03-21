//! This is a regular crate doc comment, but it also contains a partial
//! Cargo manifest.  Note the use of a *fenced* code block, and the
//! `cargo` "language".
//!
//! ```cargo
//! [dependencies]
//! csv = "1.2"
//! serde = { version = "1", features = ["derive"] }
//! serde_yaml = "0.9"
//! serde_json = "1"
//! serde_derive = "1"
//! itertools = "0.10"
//! ```

use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use csv;
use serde::{Serialize, Deserialize};
use itertools::Itertools;

#[derive(Debug, Deserialize)]
#[serde(rename_all = "SCREAMING_SNAKE_CASE")]
struct InputRecord {
  chrom: String,
  pos: u64,
  r#ref: char,
  alt: char,
  hom_ref: f64,
  het: f64,
  hom_alt: f64,
}

#[derive(Debug, Serialize)]
struct OutputRecord {
  variant_key: String,
  ml_genotype: char,
  likelihoods: String,
}

fn main() -> Result<(), Box<dyn Error>> {
  let _stdout_redirect = snakemake.redirect_stdout(&snakemake.log[0])?;
  let _stderr_redirect = snakemake.redirect_stderr(&snakemake.log[0])?;

  let mut writer = csv::WriterBuilder::new()
      .delimiter(b'\t')
      .from_path(&snakemake.output.ml)?;

  let genotype_mapping_file = File::open(&snakemake.input.genotype_mapping)?;
  let genotype_mapping: HashMap<String, char> = serde_yaml::from_reader(genotype_mapping_file)?;
  let genotype_order_file = File::open(&snakemake.input.genotype_order)?;
  let genotype_order: HashMap<char, usize> = serde_yaml::from_reader(genotype_order_file)?;

  let mut rdr = csv::ReaderBuilder::new()
      .has_headers(true)
      .delimiter(b'\t')
      .from_path(&snakemake.input.gt_likelihoods)?;
  
  fn iupac_from_alleles(allele_one: &char, allele_two: &char, lookup: &HashMap<String, char>) -> char {
    let chars = format!("{}{}", allele_one, allele_two);
    *lookup.get(&chars).unwrap()
  }

//  let mut record_in_index: HashMap<String, usize> = HashMap::new();
//  let mut record_fields_out: usize;
//  let mut record_fields_in: usize;
//  {
//      let header_in = rdr.headers()?;
//      // columns that are added for each line below
//      record_fields_in = header_in.len();
//
//      record_in_index = header_in
//          .iter()
//          .enumerate()
//          .map(|(i, name)| (name.to_string(), i))
//          .collect();
//
//  }
//
//  // rough calculation:
//  // CHROM:          15
//  // POS:             9
//  // REF:             1
//  // ALT:             1
//  // PROBs: 21 x 3 = 63
//  // ==================
//  //                 89
//  // 128 gives ample extra space, in case e.g. CHROM is much longer
//  let input_buffer_size: usize = 128;
//  let mut record_in = csv::StringRecord::with_capacity(input_buffer_size, record_fields_in);
//
//  let header_out = vec!["var_key", "ml_genotype", "likelihoods"];
//  let record_out_index: HashMap<String, usize> = header_out
//      .iter()
//      .enumerate()
//      .map(|(i, name)| (name.to_string(), i))
//      .collect();
//  writer.write_record(header_out)?;
//  // rough calculation:
//  // var_key: 15 + 9 + 1 + 1 + 3 =  29
//  // ml_genotype:                =   2
//  // PROBs:              21 x 10 = 210
//  // =================================
//  //                               241
//  // 256 gives some extra space, in case e.g. CHROM is much longer
//  let mut var_key = String::with_capacity(30);
//  let output_buffer_size = 241;
//  let mut record_out = csv::StringRecord::with_capacity(output_buffer_size, header_out.len());
 // rough calculations:
 // CHROM:          16
 // ==================
 //                 89
  let mut record_in = InputRecord {
    chrom: String::with_capacity(16),
    pos: 0,
    r#ref: 'N',
    alt: 'N',
    hom_ref: 1.0,
    het: 1.0,
    hom_alt: 1.0,
  };
  // rough calculation:
  // var_key: 16 + 9 + 1 + 1 + 3 =  30
  // ml_genotype:                =   2
  // PROBs:              21 x 10 = 210
  // =================================
  //                               241
  // 256 gives some extra space, in case e.g. CHROM is much longer
  let mut variant_key = String::with_capacity(30);
  let mut ml_genotype: char;
  let mut likelihoods_string = String::with_capacity(210);
  let mut likelihoods_vec: Vec<f64>;

  for result in rdr.deserialize() {
    // get existing entries
    record_in = result?;

    // create one single string as more efficient merge key across cells downstream
    variant_key = format!(
      "{}:{}_{}_{}",
      record_in.chrom,
      record_in.pos,
      record_in.r#ref,
      record_in.alt,
    );

    let iupacs: HashMap<String, char> = HashMap::from([
      (
        String::from("hom_ref"),
        iupac_from_alleles(&record_in.r#ref, &record_in.r#ref, &genotype_mapping)
      ),
      (
        String::from("het"),
        iupac_from_alleles(&record_in.r#ref, &record_in.alt, &genotype_mapping)
      ),
      (
        String::from("hom_alt"),
        iupac_from_alleles(&record_in.alt, &record_in.alt, &genotype_mapping)
      ),
    ]);
    
    ml_genotype = *iupacs.get(&snakemake.wildcards.genotype).unwrap();

    // reset likelihoods
    likelihoods_vec = vec![0.0; 10];
    // insert ProSolo genotype probabilities
    likelihoods_vec[*genotype_order.get(iupacs.get("hom_ref").unwrap()).unwrap()] = record_in.hom_ref;
    likelihoods_vec[*genotype_order.get(iupacs.get("het").unwrap()).unwrap()] = record_in.het;
    likelihoods_vec[*genotype_order.get(iupacs.get("hom_alt").unwrap()).unwrap()] = record_in.hom_alt;
    likelihoods_string = format!("{:.17}", likelihoods_vec.iter().format(","));

    writer.serialize(OutputRecord{
      variant_key: variant_key,
      ml_genotype: ml_genotype,
      likelihoods: likelihoods_string,
    })?;
  }
  writer.flush()?;
  Ok(())
}