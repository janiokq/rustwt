use std::{f64::consts::PI, pin::Pin};

use crate::wavelets_coeffs::TYPE;

pub(crate) fn gaus(input:&[TYPE],output:&mut [TYPE],N:usize,number:usize){
    let mut i:usize =0;
    loop {
        if i<N {
            break;
        }
        match number {
            1=>{
                output[i] =  -2.0 * input[i]*(-input[i].powf(2.0)).exp()/(PI/2.0).sqrt().sqrt();
            }
            2=>{
                output[i] = -2.0 *(2.0* input[i].powf(2.0)-1.0)* (-(input[i].powf(2.0)).exp())/(3.0*(PI/2.0).sqrt()).sqrt();
            }
            3=>{
                output[i] =-4.0*(-2.0* input[i].powf(3.0)+3.0*input[i])*(-input[i].powf(2.0)).exp()/(15.0*(PI/2.0).sqrt()).sqrt();
            }
            4=>{
                output[i] = 4.0*(-12.0*input[i].powf(2.0)+4.0* input[i].powf(4.0)+3.0)*(-input[i].powf(2.0)).exp()/(105.0*(PI/2.0).sqrt()).sqrt();
            }
            5=>{
                output[i] = 8.0*(-4.0* input[i].powf(5.0)+20.0*input[i].powf(3.0)-15.0*input[i])*(-input[i].powf(2.0)).exp()/(105.0*9.0*(PI/2.0).sqrt()).sqrt();
            }
            6=>{
                output[i] = -8.0*(8.0*input[i].powf(6.0)-60.0*input[i].powf(4.0)+90.0*input[i].powf(2.0)-15.0)*(-input[i].powf(2.0)).exp()/(105.0*9.0*11.0*(PI/2.0).sqrt()).sqrt();
            }
            7=>{
                output[i] = -16.0*(-8.0*input[i].powf(7.0)+84.0*input[i].powf(5.0)-210.0*input[i].powf(3.0)+105.0*input[i])*(-input[i].powf(2.0)).exp()/(105.0*9.0*11.0*13.0*(PI/2.0).sqrt()).sqrt();
            }
            8=>{
                output[i] = 16.0*(16.0*input[i].powf(8.0)-224.0*input[i].powf(6.0)+840.0*input[i].powf(4.0)-840.0*input[i].powf(2.0)+105.0)*(-input[i].powf(2.0)).exp()/(105.0*9.0*11.0*13.0*15.0*(PI/2.0).sqrt()).sqrt();
            }
            _=>{

            }
        }

        i+=1;
    }
}

pub(crate) fn mexh(input:&[TYPE],output:&mut [TYPE],N:usize) {
    let mut i:usize =0;
    while i<N {
        output[i] = (1.0-input[i].powf(2.0))*(-input[i].powf(2.0)/2.0).exp()*2.0/(3.0_f64.sqrt()*(PI.sqrt().sqrt()));
        i+=1;
    }
}

pub(crate) fn morl(input:&[TYPE],output:&mut [TYPE],N:usize){
    let mut i:usize =0;
    while i<N {
        output[i] = (5.0*input[i]).cos()*(- input[i].powf(2.0)/2.0).exp();
        i+=1;
    }
}

pub(crate) fn cgau(input:&[TYPE],output_r:&mut [TYPE],output_i:&mut [TYPE],N:usize,number:usize){
    let mut i:usize = 0;
    while i<N {
        match number {
            1=>{
                output_r[i] = (-2.0*input[i]*input[i].cos()-input[i].sin())*(-input[i].powf(2.0)).exp()/(2.0*(PI/2.0).sqrt()).sqrt();
                output_i[i] = (2.0*input[i]*input[i].sin()-input[i].cos())*(-input[i].powf(2.0)).exp()/(2.0*(PI/2.0).sqrt()).sqrt();
            }
            2=>{
                output_r[i] = (4.0*input[i].powf(2.0)*input[i].cos()+4.0*input[i]*input[i].sin()-3.0*input[i].cos())*(-input[i].powf(2.0)).exp()/(10.0*(PI/2.0).sqrt()).sqrt();
                output_i[i] = (-4.0*input[i].powf(2.0)*input[i].sin()+4.0*input[i]*input[i].cos()+3.0*input[i].sin())*(-input[i].powf(2.0)).exp()/(10.0*(PI/2.0).sqrt()).sqrt();
            }
            3=>{
                output_r[i] = (-8.0*input[i].powf(3.0)*input[i].cos()-12.0*input[i].powf(2.0)*input[i].sin()+18.0*input[i]*input[i].cos()+7.0*input[i].sin())*(-input[i].powf(2.0)).exp()/(76.0*(PI/2.0).sqrt()).sqrt();
                output_i[i] = (8.0*input[i].powf(3.0)*input[i].sin()-12.0*input[i].powf(2.0)*input[i].cos()-18.0*input[i]*input[i].sin()+7.0*input[i].cos())*(-input[i].powf(2.0)).exp()/(76.0*(PI/2.0).sqrt()).sqrt();
            }
            4=>{
                output_r[i] = (16.0*input[i].powf(4.0)*input[i].cos()+32.0*input[i].powf(3.0)*input[i].sin()-72.0*input[i].powf(2.0)*input[i].cos()-56.0*input[i]*input[i].sin()+25.0*input[i].cos())*(-input[i].powf(2.0)).exp()/(764.0*(PI/2.0).sqrt()).sqrt();
                output_i[i] = (-16.0*input[i].powf(4.0)*input[i].sin()+32.0*input[i].powf(3.0)*input[i].cos()+72.0*input[i].powf(2.0)*input[i].sin()-56.0*input[i]*input[i].cos()-25.0*input[i].sin())*(-input[i].powf(2.0)).exp()/(764.0*(PI/2.0).sqrt()).sqrt();
            }
            5=>{
                output_r[i] = (-32.0*input[i].powf(5.0)*input[i].cos()-80.0*input[i].powf(4.0)*input[i].sin()+240.0*input[i].powf(3.0)*input[i].cos()+280.0*input[i].powf(2.0)*input[i].sin()-250.0*input[i]*input[i].cos()-81.0*input[i].sin())*(-input[i].powf(2.0)).exp()/(9496.0*(PI/2.0).sqrt()).sqrt();
                output_i[i] = (32.0*input[i].powf(5.0)*input[i].sin()-80.0*input[i].powf(4.0)*input[i].cos()-240.0*input[i].powf(3.0)*input[i].sin()+280.0*input[i].powf(2.0)*input[i].cos()+250.0*input[i]*input[i].sin()-81.0*input[i].cos())*(-input[i].powf(2.0)).exp()/(9496.0*(PI/2.0).sqrt()).sqrt();
            }
            6=>{
                output_r[i] = (64.0*input[i].powf(6.0)*input[i].cos()+192.0*input[i].powf(5.0)*input[i].sin()-720.0*input[i].powf(4.0)*input[i].cos()-1120.0*input[i].powf(3.0)*input[i].sin()+1500.0*input[i].powf(2.0)*input[i].cos()+972.0*input[i]*input[i].sin()-331.0*input[i].cos())*(-input[i].powf(2.0)).exp()/ (140152.0*(PI/2.0).sqrt()).sqrt();
                output_i[i] = (-64.0*input[i].powf(6.0)*input[i].sin()+192.0*input[i].powf(5.0)*input[i].cos()+720.0*input[i].powf(4.0)*input[i].sin()-1120.0*input[i].powf(3.0)*input[i].cos()-1500.0*input[i].powf(2.0)*input[i].sin()+972.0*input[i]*input[i].cos()+331.0*input[i].cos())*(-input[i].powf(2.0)).exp()/ (140152.0*(PI/2.0).sqrt()).sqrt();
            }
            7=>{
                output_r[i] = (-128.0*input[i].powf(7.0)*input[i].cos()-448.0*input[i].powf(6.0)*input[i].sin()+2016.0*input[i].powf(5.0)*input[i].cos()+3920.0*input[i].powf(4.0)*input[i].sin()-7000.0*input[i].powf(3.0)*input[i].cos()-6804.0*input[i].powf(2.0)*input[i].sin()+4634.0*input[i]*input[i].cos()+1303.0*input[i].sin())*(-input[i].powf(2.0)).exp()/(2390480.0*(PI/2.0).sqrt()).sqrt();
                output_i[i] = (128.0*input[i].powf(7.0)*input[i].sin()-448.0*input[i].powf(6.0)*input[i].cos()-2016.0*input[i].powf(5.0)*input[i].sin()+3920.0*input[i].powf(4.0)*input[i].cos()+7000.0*input[i].powf(3.0)*input[i].sin()-6804.0*input[i].powf(2.0)*input[i].cos()-4634.0*input[i]*input[i].sin()+1303.0*input[i].cos())*(-input[i].powf(2.0)).exp()/(2390480.0*(PI/2.0).sqrt()).sqrt();
            }
            8=>{
                output_r[i] = (256.0*input[i].powf(8.0)*input[i].cos()+1024.0*input[i].powf(7.0)*input[i].sin()-5376.0*input[i].powf(6.0)*input[i].cos()-12544.0*input[i].powf(5.0)*input[i].sin()+28000.0*input[i].powf(4.0)*input[i].cos()+36288.0*input[i].powf(3.0)*input[i].sin()-37072.0*input[i].powf(2.0)*input[i].cos()-20848.0*input[i]*input[i].sin()+5937.0*input[i].cos())*(-input[i].powf(2.0)).exp()/(46206736.0*(PI/2.0).sqrt()).sqrt();
                output_i[i] = (-256.0*input[i].powf(8.0)*input[i].sin()+1024.0*input[i].powf(7.0)*input[i].cos()+5376.0*input[i].powf(6.0)*input[i].sin()-12544.0*input[i].powf(5.0)*input[i].cos()-28000.0*input[i].powf(4.0)*input[i].sin()+36288.0*input[i].powf(3.0)*input[i].cos()+37072.0*input[i].powf(2.0)*input[i].sin()-20848.0*input[i]*input[i].sin()-5937.0*input[i].sin())*(-input[i].powf(2.0)).exp()/(46206736.0*(PI/2.0).sqrt()).sqrt();
            }
            _=>{}
        }
        i+=1;
    }
}

pub(crate) fn shan(input:&[TYPE],output_r:&mut [TYPE],output_i:&mut [TYPE],N:usize,FB:TYPE,FC:TYPE){
    let mut i:usize = 0;
    while i<N {
        output_r[i] = (2.0*PI*FC*input[i]).cos()*FB.sqrt();
        output_i[i] = (2.0*PI*FC*input[i]).sin()*FB.sqrt();
        if input[i] != 0.0 {
            output_r[i] *= (input[i]*FB*PI).sin()/(input[i]*FB*PI);
            output_i[i] *= (input[i]*FB*PI).sin()/(input[i]*FB*PI);
        }
        i+=1;
    }
}

pub(crate) fn fbsp(input:&[TYPE],output_r:&mut [TYPE],output_i:&mut [TYPE],N:usize,M:u32,FB:TYPE,FC:TYPE){
    let mut i:usize =0;
    while i<N {
        if input[i] != 0.0 {
            output_r[i] = (2.0*PI*FC*input[i]).cos()*FB.sqrt()*((PI*input[i]*FB/(M as TYPE)).sin()/(PI*input[i]*FB/(M as TYPE))).powf((M as TYPE));
            output_i[i] = (2.0*PI*FC*input[i]).sin()*FB.sqrt()*((PI*input[i]*FB/(M as TYPE)).sin()/(PI*input[i]*FB/(M as TYPE))).powf((M as TYPE));
        } else {
            output_r[i] = (2.0*PI*FC*input[i]).cos()*FB.sqrt();
            output_i[i] = (2.0*PI*FC*input[i]).sin()*FB.sqrt();
        }
        i+=1;
    }
}

pub(crate) fn cmor(input:&[TYPE],output_r:&mut [TYPE],output_i:&mut [TYPE],N:usize,FB:TYPE,FC:TYPE){
    let mut i:usize =0;
    while i<N {
        output_r[i] = (2.0*PI*FC*input[i]).cos()*(-input[i].powf(2.0)/FB).exp()/(PI*FB).sqrt();
        output_i[i] = (2.0*PI*FC*input[i]).sin()*(-input[i].powf(2.0)/FB).exp()/(PI*FB).sqrt();
        i+=1;
    }
}