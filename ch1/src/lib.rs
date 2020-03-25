#![allow(dead_code)]


pub mod vec;
pub mod color;
pub mod canvas;
pub mod transform;
pub mod shape;
pub mod material;
pub mod camera;

pub fn floats_equal(f1:f64, f2:f64) -> bool {
    return (f1 - f2).abs() < 0.0001;
}

#[cfg(test)]
mod tests {



}
