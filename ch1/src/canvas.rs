extern crate image;


use std::io::Write;
use std::path::Path;
use std::fs::File;
use std::error::Error;
use std::convert::TryInto;
use std::time::Duration;

use super::color::Color;

use image::{Rgba, RgbaImage, Pixel, Frame, Delay};
use image::gif::Encoder;


#[derive(Debug)]
pub struct Canvas {
    width: u32,
    height: u32,
    gamma: f32,
    pixels: Vec<Color>
}

impl Canvas {

    pub fn new(w:u32, h:u32) -> Canvas {
        let capacity:usize = (w*h).try_into().unwrap();
        let mut pixels = Vec::with_capacity(capacity);
        Canvas::init_pixels(&mut pixels, w, h);
        Canvas { width:w, height:h, gamma: 1.0, pixels: pixels }
    }

    pub fn set_gamma(&mut self, gamma:f32) {
        self.gamma = gamma;
    }

    fn init_pixels(_pixels: &mut Vec<Color>, w:u32, h:u32) {
        let c = Color::new(0.0, 0.0, 0.0);
        let capacity:usize = (w*h).try_into().unwrap();
        _pixels.resize(capacity, c); 
    }

    pub fn write_pixel(&mut self, x:u32, y:u32, c:Color) {
        let idx:usize = (y*self.width+x).try_into().unwrap();
        self.pixels[idx] = c;
    }

    pub fn get_pixel(&self, x:u32, y:u32) -> Color {
        let idx:usize = (y*self.width+x).try_into().unwrap();
        self.pixels[idx]
    }

    pub fn frame_to_file<W: std::io::Write>(&self, encoder: &mut Encoder<W>) {
        let d = Duration::from_millis(75);
        encoder.encode_frame(Frame::from_parts
                             (self.to_imgbuf::<Rgba<u8>, Vec<u8>>(),
                             0, 0, Delay::from_saturating_duration(d))
        ).unwrap();
    }

    fn to_imgbuf<P: Pixel, C>(&self) -> RgbaImage {

        let mut imgbuf: RgbaImage =  RgbaImage::new(self.width, self.height);
        
        for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
            let c = self.get_pixel(x, y);
            // the casts below are needed to avoid obscure compiler errors, apparently because
            // the subpixel type (u8 here) needs to satisfy an internal trait of the image
            // library

            *pixel = image::Rgba([c.red_scaled_gamma(255, self.gamma) as u8, 
                                  c.green_scaled_gamma(255, self.gamma) as u8,
                                  c.blue_scaled_gamma(255, self.gamma) as u8,
                                  std::u8::MAX
            ]);
        }
        return imgbuf;
    }

    pub fn write_to_file(&self, file_name:&str) {

        self.to_imgbuf::<Rgba<u8>, Vec<u8>>().save(Path::new(file_name)).unwrap();
                    
    }

    pub fn write_to_file_simple(&self, file_name:&str) {
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

    pub fn dump(&self) {
        for i in 0..self.pixels.len() {
            println!("pixel at {}: {:?}", i, self.pixels[i]);
        }
    }
}

#[cfg(test)]
mod tests {

    use crate::color::*;
    use super::Canvas;

//    #[test]
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


