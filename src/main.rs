// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

mod application;

fn main() {
    std::process::exit(match application::run() {
        Ok(_) => 0,
        Err(_) => 1,
    });
}
