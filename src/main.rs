use rhilbert::*;

fn main() {
    let m = 3u8;
    let n = 2u8;
    let mut hc = HilbertCurve::new(m, n);
    let max_value = hc.max_value();
    let max_coord = hc.max_coord();
    println!("m = {}, n = {}, max_value = {}, max_coord = {}", m, n, max_value, max_coord);
    for i in 0..max_value {
        let a = hc.value_to_point(i);
        let h = hc.point_to_value(a.clone());
        println!("i = {}, a = {:?}, h = {}", i, a, h);
    }
}