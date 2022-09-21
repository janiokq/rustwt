pub enum MODE {
    MODE_INVALID = -1,
    MODE_ZEROPAD = 0, /* default, signal extended with zeros */
    MODE_SYMMETRIC,   /* signal extended symmertically (mirror)
                       * also known as half-sample symmetric
                       * For extensions greater than signal length,
                       * mirror back and forth:
                       * 2 3 3 2 1 | 1 2 3 | 3 2 1 1 2
                       */
    MODE_CONSTANT_EDGE, /* signal extended with the border value */
    MODE_SSOMOOTH,      /* linear extrapolation (first derivative) */
    MODE_PERIODIC,      /* signal is treated as being periodic */
    MODE_PERIODIZATION, /* signal is treated as being periodic, minimal output length   */
    MODE_REFLECT,       /* signal extended symmetrically (reflect)
                         * also known as whole-sample symmetric
                         * For  greater than signal length,
                         * reflect back and forth without repeating edge values:
                         * 1 2 3 2 | 1 2 3 | 2 1 2 3
                         */
    MODE_ANTISYMMMETRIC, /* antisymmetric version of "MODE_SYMMETRIC"
                          * also known as half-sample antisymmetric
                          * 2 3 -3 -2 -1 | 1 2 3 | -3 -2 -1 1 2
                          */
    MODE_ANTIREFLECT, /* antisymmetric version of "MODE_REFLECT"
                       * also know whole-sample antisymmetric
                       * 0 -1 -2 -1 0 | 1 2 3 | 4 5 6 5 4
                       */
    MODE_MAX,
}

pub fn dwt_buffer_length(input_len: usize, filter_len: usize, mode: MODE) -> usize {
    if input_len < 1 || filter_len < 1 {
        return 0;
    }
    match mode {
        MODE::MODE_PERIODIZATION => input_len / 2 + (if (input_len % 2) == 1 { 1 } else { 0 }),
        _ => (input_len + filter_len - 1) / 2,
    }
}

pub fn reconstruction_buffer_length(coeffs_len: usize, filter_len: usize) -> usize {
    if coeffs_len < 1 || filter_len < 1 {
        return 0;
    }
    2 * coeffs_len + filter_len - 2
}

pub fn idwt_buffer_length(coeffs_len: usize, filter_len: usize, mode: MODE) -> usize {
    match mode {
        MODE::MODE_PERIODIZATION => 2 * coeffs_len,
        _ => 2 * coeffs_len - filter_len + 2,
    }
}

pub fn dwt_max_level(input_len: usize, filter_len: usize) -> usize {
    if filter_len <= 1 || input_len < (filter_len - 1) {
        return 0;
    }
    f64::log2((input_len / (filter_len - 1)) as f64) as usize
}

/* check how many times input_len is divisible by 2 */
pub fn swt_max_level(input_len: usize) -> usize {
    let mut j: usize = 0;
    let mut input: usize = input_len;
    while input > 0 {
        if (input % 2) == 1 {
            return j;
        }
        input /= 2;
        j += 1;
    }
    j
}
