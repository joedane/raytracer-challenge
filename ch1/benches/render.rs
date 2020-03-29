
use criterion::{criterion_group, criterion_main, Criterion, Benchmark};

use ch1::shape::{World, Sphere, Plane};
use ch1::material::{Material, CheckerPattern};
use ch1::transform::Matrix;
use ch1::vec::{Point, Vector};
use ch1::camera::Camera;
use ch1::color::Color;

fn setup_world() -> World {
    let mut world = World::new(Default::default());
    
    let shape = Sphere::new_with_transform_and_material(1.0, 
                                                       Matrix::identity().translation(-0.5, 1., 0.5),
                                                        Material { color:Some(Color::RED),
                                                                   diffuse:0.1,
                                                                   transparency:1.0,
                                                                   refractive_index:1.15,
                                                                   specular:0.1,
                                                                   ambient:0.1,
                                                                  ..Material::DEFAULT
                                                       });
    world.add_shape(Box::new(shape));

    let shape = Sphere::new_with_transform_and_material(1.0,
                                                        Matrix::identity().scaling(0.5, 0.5, 0.5).translation(1., 0.7, -3.5),
                                                        Material { color:Some(Color::GREEN),
                                                                   diffuse:0.1,
                                                                   transparency:1.0,
                                                                   ambient:0.1,
                                                                   refractive_index:1.5,
                                                                   specular:0.1,
                                                                   ..Material::DEFAULT
                                                        });
    
    world.add_shape(Box::new(shape));

    let shape = Sphere::new_with_transform_and_material(1.0,
                                                        Matrix::identity()
                                                        .scaling(0.8, 0.8, 0.8)
                                                        .translation(-2.5, 0.53, -0.75).
                                                        rotation_x(std::f64::consts::PI/4.0),
                                                        Material { color:Some(Color::BLUE),
                                                                   diffuse:0.7,
                                                                   specular:0.3,
                                                                   ..Material::DEFAULT
                                                        });
    
    world.add_shape(Box::new(shape));

    let plane = Plane::new_with_transform_and_material(Matrix::identity().translation(0., -3., 0.), 
                                                       Material { pattern:Some(Box::new(CheckerPattern::new(Color::WHITE, Color::BLACK))),
                                                                  diffuse:0.2,
                                                                  ambient:0.6,
                                                                  specular:0.3,
                                                                  ..Material::DEFAULT
                                                        }); 
    world.add_shape(Box::new(plane));

    return world;

}

fn render() {
    let camera = Camera::new_with_transform(400, 300, std::f64::consts::PI/3.0,
                                            Matrix::make_view_transform(Point::new(0., 7., 0.),
                                                                        Point::new(0., 0., 0.),
                                                                        Vector::new(1., 0., 0.)));
    let world = setup_world();
    let canvas = camera.render(&world);
}

fn render_async() {
    let camera = Camera::new_with_transform(400, 300, std::f64::consts::PI/3.0,
                                            Matrix::make_view_transform(Point::new(0., 7., 0.),
                                                                        Point::new(0., 0., 0.),
                                                                        Vector::new(1., 0., 0.)));
    let world = setup_world();
    let canvas = camera.render_async(&world);
}

fn criterion_benchmark(c: &mut Criterion) {
    c.bench("simple render",
            Benchmark::new("f1", |b| b.iter(|| render_async())).sample_size(10)
    );
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
