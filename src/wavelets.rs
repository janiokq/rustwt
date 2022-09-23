use std::{mem::swap};

use crate::wavelets_coeffs::{coif, db, sqrt2, sym, bior, dmey};

macro_rules! nelems {
    ($x:expr) => {
        std::mem::size_of_val(&$x) / std::mem::size_of_val(&$x[0])
    };
}

/* Wavelet symmetry properties */
pub(crate) enum SYSMMETRY {
    UNKNOWN = -1,
    ASYMMETRIC = 0,
    NEAR_SYMMETRIC = 1,
    SYMMETRIC = 2,
    ANTI_SYMMETRIC = 3,
}

/* Wavelet name */

pub(crate) enum WAVELET_NAME {
    HAAR,
    RBIO,
    DB,
    SYM,
    COIF,
    BIOR,
    DMEY,
    GAUS,
    MEXH,
    MORL,
    CGAU,
    SHAN,
    FBSP,
    CMOR,
}

/* Wavelet structure holding pointers to filter arrays and property attributes */
pub(crate) struct BaseWavelet {
    /* Wavelet properties */
    support_width: i32,
    symmetry: SYSMMETRY,
    orthogonal: u32,
    biorthogonal: u32,
    compact_support: u32,

    _builtin: i32,
    family_name: String,
    short_name: String,
}


pub(crate) struct DiscreteWavelet {
    base: BaseWavelet,
    dec_hi_double: f64, /* highpass decomposition */
    dec_lo_double: f64, /* lowpass decomposition */
    rec_hi_double: f64, /* highpass reconstruction */
    rec_lo_double: f64, /* lowpass reconstruction */

    pub dec_len: usize, /* length of decomposition filter */
    pub rec_len: usize, /* length of reconstruction filter */

    vanishing_moments_psi: i32,
    vanishing_moments_phi: i32,
}


pub(crate) struct ContinuousWavelet {
    base: BaseWavelet,
    lower_bound: f32,
    upper_bound: f32,
    /* Parameters for shan, fbsp, cmor */
    complex_cwt: i32,
    center_frequency: f32,
    bandwidth_frequency: f32,
    fbsp_order: u32,
}
impl  ContinuousWavelet {
    pub fn new() -> Self{
        ContinuousWavelet{
            base: BaseWavelet::new(),
            lower_bound: Default::default(),
            upper_bound: Default::default(),
            complex_cwt: Default::default(),
            center_frequency: Default::default(),
            bandwidth_frequency: Default::default(),
            fbsp_order: Default::default(),
        }
    }
}

pub(crate) fn is_discrete_wavelet(name: WAVELET_NAME) -> i32 {
    match name {
        WAVELET_NAME::HAAR => return 0,
        WAVELET_NAME::RBIO => return 1,
        WAVELET_NAME::DB => return 1,
        WAVELET_NAME::SYM => return 1,
        WAVELET_NAME::COIF => return 1,
        WAVELET_NAME::BIOR => return 1,
        WAVELET_NAME::DMEY => return 1,
        WAVELET_NAME::GAUS => return 0,
        WAVELET_NAME::MEXH => return 0,
        WAVELET_NAME::MORL => return 0,
        WAVELET_NAME::CGAU => return 0,
        WAVELET_NAME::SHAN => return 0,
        WAVELET_NAME::FBSP => return 0,
        WAVELET_NAME::CMOR => return 0,
        _ => return 0,
    }
}

impl BaseWavelet {
    pub fn new() -> Self {
        BaseWavelet {
            support_width: -1,
            symmetry: SYSMMETRY::UNKNOWN,
            orthogonal: Default::default(),
            biorthogonal: Default::default(),
            compact_support: Default::default(),
            _builtin: Default::default(),
            family_name: Default::default(),
            short_name: Default::default(),
        }
    }
}

impl DiscreteWavelet {
    pub fn new() -> Self {
        DiscreteWavelet {
            base: BaseWavelet::new(),
            dec_hi_double: Default::default(),
            dec_lo_double: Default::default(),
            rec_hi_double: Default::default(),
            rec_lo_double: Default::default(),
            dec_len: Default::default(),
            rec_len: Default::default(),
            vanishing_moments_psi: Default::default(),
            vanishing_moments_phi: Default::default(),
        }
    }
}

/*
 * Allocate Wavelet struct and set its attributes
 * name - (currently) a character codename of a wavelet family
 * order - order of the wavelet (ie. coif3 has order 3)
 */
pub(crate) fn discrete_wavelet(name: WAVELET_NAME, order: usize) -> Option<DiscreteWavelet> {
    match name {
        /* Haar wavelet */
        WAVELET_NAME::HAAR => {
            /* the same as db1 */
            let mut w = discrete_wavelet(WAVELET_NAME::DB, 1);
            if let Some(mut wav) = w {
                wav.base.family_name = "Haar".to_string();
                wav.base.short_name = "haar".to_string();
                Some(wav)
            } else {
                None
            }
        }
        /* Reverse biorthogonal wavelets family */
        WAVELET_NAME::RBIO => {
            /* rbio is like bior, only with switched filters */
            let mut w = discrete_wavelet(WAVELET_NAME::BIOR, order);
            if let Some(mut wav) = w {
                swap(&mut wav.dec_len, &mut wav.rec_len);
                swap(&mut wav.rec_lo_double, &mut wav.dec_lo_double);
                swap(&mut wav.rec_hi_double, &mut wav.dec_hi_double);
                // {
                //     let mut i = 0;
                //     let mut j = wav.rec_len -1;
                //     loop {
                //         if !(i<j) {
                //             break;
                //         }
                //         i+=1;
                //         j-=1;
                //     }
                // }
                wav.base.family_name = "Reverse biorthogonal".to_string();
                wav.base.short_name = "rbio".to_string();
                Some(wav)
            } else {
                None
            }
        }
        /* Daubechies wavelets family */
        WAVELET_NAME::DB => {
            let coeff_idx = order - 1;
            if coeff_idx >= nelems!(db) {
                return None;
            }
            let mut w = blank_discrete_wavelet(2 * order);
            if let Some(mut wav) = w {
                wav.vanishing_moments_psi = order as i32;
                wav.vanishing_moments_phi = 0;
                wav.base.support_width = (2 * order) as i32;
                wav.base.orthogonal = 1;
                wav.base.biorthogonal = 1;
                wav.base.symmetry = SYSMMETRY::ASYMMETRIC;
                wav.base.compact_support = 1;
                wav.base.family_name = "Daubechies".to_string();
                wav.base.short_name = "db".to_string();
                {
                    let mut i = 0;
                    loop {
                        if !(i < wav.rec_len) {
                            break;
                        }
                        wav.rec_lo_double = db[coeff_idx][i];
                        wav.dec_lo_double = db[coeff_idx][wav.dec_len - 1 - i];
                        let mut rec_hi_double: f64 = 0.0;
                        if (i % 2) == 1 {
                            rec_hi_double = -1.0;
                        } else {
                            rec_hi_double = 1.0;
                        }
                        wav.rec_hi_double = rec_hi_double * (db[coeff_idx][wav.dec_len - 1 - i]);
                        let mut dec_hi_double: f64 = 0.0;
                        if ((wav.dec_len - 1 - i) % 2) == 1 {
                            dec_hi_double = -1.0;
                        } else {
                            dec_hi_double = 1.0;
                        }
                        wav.dec_hi_double = dec_hi_double * db[coeff_idx][i];

                        i += 1;
                    }
                }
                None
            } else {
                None
            }
        }
        /* Symlets wavelets family */
        WAVELET_NAME::SYM => {
            let coeff_idx = order - 2;
            if coeff_idx >= nelems!(sym) {
                return None;
            }
            let mut w = blank_discrete_wavelet(2 * order);
            if let Some(mut wav) = w {
                wav.vanishing_moments_psi = order as i32;
                wav.vanishing_moments_phi = 0;
                wav.base.support_width = (2 * order - 1) as i32;
                wav.base.orthogonal = 1;
                wav.base.biorthogonal = 1;
                wav.base.symmetry = SYSMMETRY::NEAR_SYMMETRIC;
                wav.base.compact_support = 1;
                wav.base.family_name = "Symlets".to_string();
                wav.base.short_name = "sym".to_string();
                {
                    let mut i = 0;
                    loop {
                        if !(i < wav.rec_len) {
                            break;
                        }
                        wav.rec_lo_double = sym[coeff_idx][i];
                        wav.dec_lo_double = sym[coeff_idx][wav.dec_len - 1 - i];
                        let mut rec_hi_double: f64 = 0.0;
                        if (i % 2) == 1 {
                            rec_hi_double = -1.0;
                        } else {
                            rec_hi_double = 1.0;
                        }
                        wav.rec_hi_double = rec_hi_double * (sym[coeff_idx][wav.dec_len - 1 - i]);
                        let mut dec_hi_double: f64 = 0.0;
                        if ((wav.dec_len - 1 - i) % 2) == 1 {
                            dec_hi_double = -1.0;
                        } else {
                            dec_hi_double = 1.0;
                        }
                        wav.dec_hi_double = dec_hi_double * sym[coeff_idx][i];
                        i += 1;
                    }
                }
                None
            } else {
                None
            }
        }
        /* Coiflets wavelets family */
        WAVELET_NAME::COIF => {
            let coeff_idx = order - 1;
            if coeff_idx >= nelems!(coif) {
                return None;
            }
            let mut w = blank_discrete_wavelet(6 * order);
            if let Some(mut wav) = w {
                wav.vanishing_moments_psi = (2 * order) as i32;
                wav.vanishing_moments_phi = (2 * order - 1) as i32;
                wav.base.support_width = (6 * order - 1) as i32;
                wav.base.orthogonal = 1;
                wav.base.biorthogonal = 1;
                wav.base.symmetry = SYSMMETRY::NEAR_SYMMETRIC;
                wav.base.compact_support = 1;
                wav.base.family_name = "Coiflets".to_string();
                wav.base.short_name = "coif".to_string();
                {
                    let mut i = 0;
                    loop {
                        if !(i < wav.rec_len) {
                            break;
                        }
                        wav.rec_lo_double = coif[coeff_idx][i] * sqrt2;
                        wav.dec_lo_double = coif[coeff_idx][wav.dec_len - 1 - i] * sqrt2;
                        let mut rec_hi_double: f64 = 0.0;
                        if (i % 2) == 1 {
                            rec_hi_double = -1.0;
                        } else {
                            rec_hi_double = 1.0;
                        }
                        wav.rec_hi_double =
                            rec_hi_double * (coif[coeff_idx][wav.dec_len - 1 - i]) * sqrt2;
                        let mut dec_hi_double: f64 = 0.0;
                        if ((wav.dec_len - 1 - i) % 2) == 1 {
                            dec_hi_double = -1.0;
                        } else {
                            dec_hi_double = 1.0;
                        }
                        wav.dec_hi_double = dec_hi_double * sym[coeff_idx][i] * sqrt2;
                        i += 1;
                    }
                }
                None
            } else {
                None
            }
        }
        /* Biorthogonal wavelets family */
        WAVELET_NAME::BIOR => {
            let mut N = order / 10;
            let mut M = order % 10;
            let mut M_idx: usize = 0;
            let mut M_max: usize = 0;
            match N {
                1 => {
                    if (M % 2 != 1 || M > 5) {
                        return None;
                    }
                    M_idx = M / 2;
                    M_max = 5;
                }
                2 => {
                    if (M % 2 != 0 || M < 2 || M > 8) {
                        return None;
                    }
                    M_idx = M / 2 - 1;
                    M_max = 8;
                }
                3 => {
                    if (M % 2 != 1) {
                        return None;
                    }
                    M_idx = M / 2;
                    M_max = 9;
                }
                4 => {}
                5 => {
                    if M != N {
                        return None;
                    }
                    M_idx = 0;
                    M_max = M;
                }
                6 => {
                    if M != 8 {
                        return None;
                    }
                    M_idx  =0;
                    M_max = 8;
                }
                _ => return None,
            }

            let mut filters_length = 0;
            if N ==1{
               filters_length =  2 * M;
            }else{
                filters_length = 2*M+2;
            }

            let mut w = blank_discrete_wavelet(filters_length);
            if let Some(mut wav) = w {
                wav.vanishing_moments_psi = (order/10) as i32;
                wav.vanishing_moments_phi = (order%10) as i32;
                wav.base.support_width = -1;
                wav.base.orthogonal = 0;
                wav.base.biorthogonal = 1;
                wav.base.symmetry = SYSMMETRY::SYMMETRIC;
                wav.base.compact_support = 1;
                wav.base.family_name = "Biorthogonal".to_string();
                wav.base.short_name = "bior".to_string();
                {
                    let mut i = 0;
                    let mut n = M_max-M;

                    loop {
                        if !(i < wav.rec_len) {
                            break;
                        }
                        wav.rec_lo_double = bior[N-1][0][i+n];
                        wav.dec_lo_double = bior[N-1][M_idx+1][wav.dec_len-1-i];
                        let mut rec_hi_double: f64 = 0.0;
                        if (i % 2) == 1 {
                            rec_hi_double = -1.0;
                        } else {
                            rec_hi_double = 1.0;
                        }
                        wav.rec_hi_double =
                            rec_hi_double * bior[N-1][M_idx+1][wav.dec_len-1-i];
                        let mut dec_hi_double: f64 = 0.0;
                        if ((wav.dec_len - 1 - i) % 2) == 1 {
                            dec_hi_double = -1.0;
                        } else {
                            dec_hi_double = 1.0;
                        }
                        wav.dec_hi_double = dec_hi_double * bior[N-1][0][i+n];
                        i += 1;
                    }
                }
                None
            } else {
                None
            }
        }
         /* Discrete FIR filter approximation of Meyer wavelet */
         WAVELET_NAME::DMEY => {
            let mut w = blank_discrete_wavelet(62);
            if let Some(mut wav) = w {
                wav.vanishing_moments_psi = -1;
                wav.vanishing_moments_phi = -1;
                wav.base.support_width = -1;
                wav.base.orthogonal = 1;
                wav.base.biorthogonal = 1;
                wav.base.symmetry = SYSMMETRY::SYMMETRIC;
                wav.base.compact_support = 1;
                wav.base.family_name = "Discrete Meyer (FIR Approximation)".to_string();
                wav.base.short_name = "dmey".to_string();
                {
                    let mut i = 0;
                    loop {
                        if !(i < wav.rec_len) {
                            break;
                        }
                        wav.rec_lo_double = dmey[i];
                        wav.dec_lo_double = dmey[wav.dec_len-1-i];
                        let mut rec_hi_double: f64 = 0.0;
                        if (i % 2) == 1 {
                            rec_hi_double = -1.0;
                        } else {
                            rec_hi_double = 1.0;
                        }
                        wav.rec_hi_double =
                            rec_hi_double * dmey[wav.dec_len-1-i];
                        let mut dec_hi_double: f64 = 0.0;
                        if ((wav.dec_len - 1 - i) % 2) == 1 {
                            dec_hi_double = -1.0;
                        } else {
                            dec_hi_double = 1.0;
                        }
                        wav.dec_hi_double = dec_hi_double * dmey[i];
                        i += 1;
                    }
                }
                None
            } else {
                None
            }
        }
        _ => None,
    }
}


pub(crate) fn continuous_wavelet(name:WAVELET_NAME,order: usize) -> Option<ContinuousWavelet> {

    match name {
        /* Gaussian wavelets */
        WAVELET_NAME::GAUS =>{
            if(order >8 ){return None;}
            let mut w = blank_continuous_wavelet();
            w.base.support_width = -1;
            w.base.orthogonal = 0;
            w.base.biorthogonal = 0;
            if order % 2 == 0 {
                w.base.symmetry = SYSMMETRY::SYMMETRIC;
            } else {
                w.base.symmetry = SYSMMETRY::ANTI_SYMMETRIC;
            }
            w.base.compact_support = 0;
            w.base.family_name = "Gaussian".to_string();
            w.base.short_name = "gaus".to_string();
            w.complex_cwt = 0;
            w.lower_bound = -5.0;
            w.center_frequency = 0.0;
            w.bandwidth_frequency = 0.0;
            w.fbsp_order = 0;
            Some(w)
        }
        WAVELET_NAME::MEXH => {
            let mut w = blank_continuous_wavelet();
            w.base.support_width = -1;
            w.base.orthogonal = 0;
            w.base.biorthogonal =0;
            w.base.symmetry =SYSMMETRY::SYMMETRIC;
            w.base.compact_support =0;
            w.base.family_name = "Mexican hat wavelet".to_string();
            w.base.short_name = "mexh".to_string();
            w.complex_cwt = 0;
            w.lower_bound = -8.0;
            w.upper_bound = 8.0;
            w.center_frequency =0.0;
            w.bandwidth_frequency =0.0;
            w.fbsp_order =0;
            Some(w)
        }
        WAVELET_NAME::MORL=>{
            let mut w = blank_continuous_wavelet();
            w.base.support_width  =-1;
            w.base.orthogonal =0;
            w.base.biorthogonal = 0;
            w.base.symmetry = SYSMMETRY::SYMMETRIC;
            w.base.compact_support =0;
            w.base.family_name = "Morlet wavelet".to_string();
            w.base.short_name = "morl".to_string();
            w.complex_cwt = 0;
            w.lower_bound = -8.0;
            w.upper_bound = 8.0;
            w.center_frequency = 0.0;
            w.bandwidth_frequency =0.0;
            w.fbsp_order = 0;
            Some(w)
        }
        WAVELET_NAME::CGAU=>{
            if order > 8 {return None;}
            let mut w = blank_continuous_wavelet();
            w.base.support_width = -1;
            w.base.orthogonal = 0;
            w.base.biorthogonal =0;
            if order % 2 ==0 {
                w.base.symmetry = SYSMMETRY::SYMMETRIC;
            } else {
                w.base.symmetry = SYSMMETRY::ANTI_SYMMETRIC;
            }
            w.base.compact_support = 0;
            w.base.family_name = "Complex Gaussian wavelets".to_string();
            w.base.short_name = "cgau".to_string();
            w.complex_cwt = 1;
            w.lower_bound = -5.0;
            w.upper_bound = 5.0;
            w.center_frequency = 0.0;
            w.bandwidth_frequency = 0.0;
            w.fbsp_order =0;
            Some(w)
        }
        WAVELET_NAME::FBSP=>{
            let mut w = blank_continuous_wavelet();
            w.base.support_width = -1;
            w.base.orthogonal =0;
            w.base.biorthogonal = 0;
            w.base.symmetry = SYSMMETRY::ASYMMETRIC;
            w.base.compact_support = 0;
            w.base.family_name = "Frequency B-Spline wavelets".to_string();
            w.base.short_name = "fbsp".to_string();
            w.complex_cwt =1;
            w.lower_bound = -20.0;
            w.upper_bound =20.0;
            w.center_frequency = 0.5;
            w.bandwidth_frequency = 1.0;
            w.fbsp_order = 2;
            Some(w)
        }
        WAVELET_NAME::CMOR=>{
            let mut w = blank_continuous_wavelet();
            w.base.support_width = -1;
            w.base.orthogonal = 0;
            w.base.biorthogonal =0;
            w.base.symmetry = SYSMMETRY::ASYMMETRIC;
            w.base.compact_support = 0;
            w.base.family_name = "Complex Morlet wavelets".to_string();
            w.base.short_name = "cmor".to_string();
            w.complex_cwt = 1;
            w.lower_bound = -8.0;
            w.upper_bound  = 8.0;
            w.center_frequency = 0.5;
            w.bandwidth_frequency = 1.0;
            w.fbsp_order =0;

            Some(w)
        }
        _=>None
    }
}

/*
 * Allocate blank Discrete Wavelet with zero-filled filters of given length
 */
pub(crate) fn blank_discrete_wavelet(filters_length: usize) -> Option<DiscreteWavelet> {
    /* pad to even length */
    let mut filters_length = filters_length;
    if filters_length > 0 && (filters_length % 2) == 1 {
        filters_length += 1;
    }

    let mut w = DiscreteWavelet::new();
    w.dec_len = filters_length;
    w.rec_len = filters_length;

    // if filters_length > 0 {
    // } else {
    // }

    w.base.support_width =-1;
    w.base.orthogonal =0;
    w.base.biorthogonal =0;
    w.base.symmetry = SYSMMETRY::UNKNOWN;
    w.base.compact_support = 0;
    w.base.family_name = "".to_string();
    w.base.short_name = "".to_string();
    w.vanishing_moments_phi =0;
    w.vanishing_moments_psi =0;
    Some(w)
}

pub(crate) fn blank_continuous_wavelet() -> ContinuousWavelet {
    let mut w = ContinuousWavelet::new();
    /* set properties to "blank" value */
    w.center_frequency = -1.0;
    w.bandwidth_frequency = -1.0;
    w.fbsp_order = 0;
    w
}