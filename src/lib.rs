mod common;
mod macro_wavelets;
mod wavelets;
mod wavelets_coeffs;

pub fn add(left: usize, right: usize) -> usize {
    left + right
}
const xxx: [i32; 3] = [2, 321, 31];

#[cfg(test)]
mod tests {
    use std::mem::swap;

    use crate::wavelets::{blank_discrete_wavelet, WAVELET_NAME};

    use super::*;

    #[test]
    fn it_works() -> () {
        let result = add(2, 2);
        assert_eq!(result, 4);
        let mut w = blank_discrete_wavelet(1);
        let s = "xx".to_string().clone();
        if let Some(mut wav) = w {
            wav.dec_len = 3;
            println!("{:?},{:?}", wav.dec_len, wav.rec_len);
            swap(&mut wav.dec_len, &mut wav.rec_len);
            println!("{:?},{:?}", wav.dec_len, wav.rec_len);
        } else {
            println!("没有创建成功");
        }

        {
            println!("测试哈哈哈哈😄")
        }
    }

    #[test]
    fn test_dwt_buffer_length() {
        use common::dwt_buffer_length;
        let x = common::dwt_max_level(10, 2);
        println!("value {}", x);
        let x = vec![12, 3, 321, 3];
        const _f: [&[i32]; 2] = [&[1], &[2, 3]];
    }
}
