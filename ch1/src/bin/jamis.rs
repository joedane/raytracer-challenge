
use std::path::Path;

fn main() {
    match ch1::lua::render_lua(Path::new("jamis.lua")) {
        Ok(s) => println!("{}", s),
        Err(e) => panic!(e)
    }    
}