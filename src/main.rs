fn main() {
    if let Err(err) = seqcheck::run() {
        eprintln!("error: {err:#}");
        std::process::exit(1);
    }
}
