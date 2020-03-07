

use std::io::Write;
use std::path::Path;
use std::fs::File;
use std::error::Error;
use super::color::Color;


#[derive(Debug)]
pub struct Canvas {
    width: usize,
    height: usize,
    pixels: Vec<Color>
}

impl Canvas {

    pub fn new(w:usize, h:usize) -> Canvas {
        let mut pixels = Vec::with_capacity(w*h);
        Canvas::init_pixels(&mut pixels, w, h);
        Canvas { width:w, height:h, pixels: pixels }
    }

    fn init_pixels(_pixels: &mut Vec<Color>, w:usize, h:usize) {
        let c = Color::new(0.0, 0.0, 0.0);
        _pixels.resize(w*h, c); 
    }

    pub fn write_pixel(&mut self, x:usize, y:usize, c:Color) {
        self.pixels[y*self.width+x] = c;
    }

    pub fn get_pixel(&self, x:usize, y:usize) -> Color {
        self.pixels[y*self.width+x]
    }

    pub fn write_to_file(&self, file_name:&str) {
        let mut outfile = match File::create(Path::new(file_name)) {
            Err(e) => panic!("Could not open output file '{}': {}", file_name, e.description()),

            Ok(file) => file,
        };
        let scale = 255;

        write!(outfile, "P3\n{} {}\n255\n", self.width, self.height).unwrap();

        // TODO - break long lines
        for i in 0..self.height {
            let mut points_on_line = 0;
            for j in 0..self.width {
                let c = self.get_pixel(j, i);
                if points_on_line > 0 {
                    write!(outfile, " ").unwrap();
                }
                write!(outfile, "{} {} {}", c.red_scaled(scale), c.green_scaled(scale), c.blue_scaled(scale)).unwrap();
                points_on_line += 1;
            }
            write!(outfile, "\n").unwrap();
        }
    }
}

#[cfg(test)]
mod tests {

    use crate::color::*;
    use super::Canvas;

    #[test]
    fn test_file_write() {
        let mut canvas = Canvas::new(3, 4);
        
        for i in 0..canvas.height {
            for j in 0..canvas.width {
                if (i*3+j) % 3 == 0 {
                    canvas.write_pixel(j, i, Color::WHITE);
                }
                else {
                    canvas.write_pixel(j, i, Color::BLACK);
                }
            }
        }
        
        canvas.write_to_file("test.txt");
    }
}


