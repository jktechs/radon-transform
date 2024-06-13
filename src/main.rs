use png::{BitDepth, ColorType};
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

fn main() {
    let (file, size) = read_png("./scan.png");
    let file: Vec<_> = file
        .into_iter()
        .map(|x| x as f64 / 255.)
        .map(|x| Complex(x, 0.))
        .collect();
    let fourier_1d: Vec<_> = (0..size)
        .flat_map(|i| ditfft2(&file[(i * size * 3)..], size, 3, false))
        .collect();
    write_png(
        "./1_fourier_1d.png",
        &complex_to_rgb(&fourier_1d),
        size,
        ColorType::Rgb,
    );
    let mut fourier_2d = vec![Complex(0., 0.); size * size];
    for i in 0..size {
        let tmp = ditfft2(&fourier_1d[i..], size, size, false);
        for j in 0..size {
            fourier_2d[j * size + i] = tmp[j];
        }
    }
    write_png(
        "./2_fourier_2d.png",
        &complex_to_rgb(&fourier_2d),
        size,
        ColorType::Rgb,
    );
    let mut transform = vec![Complex(0., 0.); size * size];
    for i in 0..size {
        for j in 0..size {
            let x = if (j as f64) < (size as f64 * 0.5) {
                j as f64 + (size as f64) * 0.5
            } else {
                j as f64 - (size as f64) * 0.5
            };
            let y = if (i as f64) < (size as f64 * 0.5) {
                i as f64 + (size as f64) * 0.5
            } else {
                i as f64 - (size as f64) * 0.5
            };
            transform[i * size + j] = fourier_2d[y as usize * size + x as usize];
        }
    }
    write_png(
        "./3_transform.png",
        &complex_to_rgb(&transform),
        size,
        ColorType::Rgb,
    );
    let mut transform_new = vec![Complex(0., 0.); size * size];
    for i in 0..size {
        for j in 0..size {
            let angle = std::f64::consts::PI * i as f64 / (size - 1) as f64;
            let x = (j as f64 - size as f64 * 0.5) * angle.cos() + size as f64 * 0.5;
            let y = (j as f64 - size as f64 * 0.5) * angle.sin() + size as f64 * 0.5;
            transform_new[i * size + j] = transform[y as usize * size + x as usize];
        }
    }
    write_png(
        "./4_transform_new.png",
        &complex_to_rgb(&transform_new),
        size,
        ColorType::Rgb,
    );
    let mut transform_new_new = vec![Complex(0., 0.); size * size];
    for i in 0..size {
        for j in 0..size {
            let x = if (j as f64) < (size as f64 * 0.5) {
                j as f64 + (size as f64) * 0.5
            } else {
                j as f64 - (size as f64) * 0.5
            };
            transform_new_new[i * size + j] = transform_new[i * size + x as usize];
        }
    }
    write_png(
        "./5_transform_new_new.png",
        &complex_to_rgb(&transform_new_new),
        size,
        ColorType::Rgb,
    );
    let radon: Vec<_> = (0..size)
        .flat_map(|i| ditfft2(&transform_new_new[(i * size)..], size, 1, true))
        .collect();
    write_png(
        "./6_radon.png",
        &complex_to_rgb(&radon),
        size,
        ColorType::Rgb,
    );
}
fn read_png(path: &str) -> (Vec<u8>, usize) {
    let decoder = png::Decoder::new(File::open(path).unwrap());
    let mut reader = decoder.read_info().unwrap();
    assert_eq!(
        reader.output_color_type(),
        (ColorType::Rgb, BitDepth::Eight)
    );
    // Allocate the output buffer.
    let mut buf = vec![0; reader.output_buffer_size()];
    // Read the next frame. An APNG might contain multiple frames.
    let info = reader.next_frame(&mut buf).unwrap();
    assert_eq!(info.width, info.height);
    let size = info.width as usize;
    assert_eq!(info.buffer_size(), 3 * size * size);
    // Inspect more details of the last read frame.
    assert!(reader.info().frame_control.is_none());
    assert_eq!(info.buffer_size(), buf.len());
    (buf, size)
}
fn write_png(str: &str, img: &[u8], size: usize, color: ColorType) {
    let path = Path::new(str);
    let file = File::create(path).unwrap();
    let w = BufWriter::new(file);

    let mut encoder = png::Encoder::new(w, size as u32, size as u32); // Width is 2 pixels and height is 1.
    encoder.set_color(color);
    encoder.set_depth(png::BitDepth::Eight);
    encoder.set_source_gamma(png::ScaledFloat::from_scaled(45455)); // 1.0 / 2.2, scaled by 100000
    encoder.set_source_gamma(png::ScaledFloat::new(1.0 / 2.2)); // 1.0 / 2.2, unscaled, but rounded
    let source_chromaticities = png::SourceChromaticities::new(
        // Using unscaled instantiation here
        (0.31270, 0.32900),
        (0.64000, 0.33000),
        (0.30000, 0.60000),
        (0.15000, 0.06000),
    );
    encoder.set_source_chromaticities(source_chromaticities);
    let mut writer = encoder.write_header().unwrap();

    writer.write_image_data(img).unwrap(); // Save
}
#[derive(Debug, Clone, Copy)]
struct Complex(f64, f64);

fn ditfft2(x: &[Complex], n: usize, s: usize, inv: bool) -> Vec<Complex> {
    if n == 1 {
        vec![x[0]]
    } else {
        let mut buffer = ditfft2(x, n / 2, 2 * s, inv);
        buffer.extend(ditfft2(&x[s..], n / 2, 2 * s, inv));
        for j in 0..(n / 2) {
            let lower = buffer[j];
            let upper = buffer[j + n / 2];
            let angle = std::f64::consts::PI * 2. * j as f64 / n as f64;
            let upper = Complex(
                angle.cos() * upper.0 + angle.sin() * upper.1 * if inv { -1. } else { 1. },
                angle.cos() * upper.1 - angle.sin() * upper.0 * if inv { -1. } else { 1. },
            );
            buffer[j] = Complex(lower.0 + upper.0, lower.1 + upper.1);
            buffer[j + n / 2] = Complex(lower.0 - upper.0, lower.1 - upper.1);
        }
        buffer
    }
}
fn complex_to_rgb(img: &[Complex]) -> Vec<u8> {
    let avrg = img
        .iter()
        .map(|x| (x.0 * x.0 + x.1 * x.1).sqrt())
        .sum::<f64>()
        / img.len() as f64;
    let dev = img
        .iter()
        .map(|x| (x.0 * x.0 + x.1 * x.1).sqrt())
        .map(|x| (x - avrg).powi(2))
        .sum::<f64>()
        .sqrt()
        / (img.len() as f64).sqrt();
    let max = avrg + 1.5 * dev;

    img.iter()
        .flat_map(|p| {
            fn hue_to_rgb(p: f64, q: f64, mut t: f64) -> f64 {
                if t < 0. {
                    t += 1.;
                }
                if t > 1. {
                    t -= 1.;
                }
                if t < 1. / 6. {
                    return p + (q - p) * 6. * t;
                }
                if t < 1. / 2. {
                    return q;
                }
                if t < 2. / 3. {
                    return p + (q - p) * (2. / 3. - t) * 6.;
                }
                p
            }
            let h = p.1.atan2(p.0) * 0.5 / std::f64::consts::PI;
            let l = (p.0 * p.0 + p.1 * p.1).sqrt() / max;
            let q = if l < 0.5 { l * 2. } else { 1. };
            let p = 2. * l - q;
            let r = hue_to_rgb(p, q, h + 1. / 3.);
            let g = hue_to_rgb(p, q, h);
            let b = hue_to_rgb(p, q, h - 1. / 3.);
            [
                (r.clamp(0., 1.) * 255.) as u8,
                (g.clamp(0., 1.) * 255.) as u8,
                (b.clamp(0., 1.) * 255.) as u8,
            ]
        })
        .collect()
}
