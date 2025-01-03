// Released under MIT License.
// Copyright (c) 2024-2025 Ladislav Bartos

mod application;

fn main() {
    std::process::exit(if application::run() { 0 } else { 1 });
}
