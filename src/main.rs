use rhilbert::*;

fn main() {
    let m = 3u8;
    let n = 2u8;
    let max_value = max_value(m, n);
    let max_coord = max_coord(m);
    println!("m = {}, n = {}, max_value = {}, max_coord = {}", m, n, max_value, max_coord);
    for i in 0..max_value {
        let a = hilbert_value_to_point(n, m, i);
        let h = point_to_hilbert_value(m, a.clone());
        println!("i = {}, a = {:?}, h = {}", i, a, h);
    }
}