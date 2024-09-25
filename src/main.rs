// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

use std::process;

mod application;

fn main() {
    if let Err(_) = application::run() {
        process::exit(1);
    }

    process::exit(0);
}
