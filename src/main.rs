// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

use std::process;

mod application;

fn main() {
    if let Err(e) = application::run() {
        eprintln!("{}", e);
        process::exit(1);
    }

    process::exit(0);
}
