

use std::char::ToUppercase;

use crate::{wavelets_coeffs::TYPE, common::MODE};

/* This file contains several functions for computing the convolution of a
 * signal with a filter. The general scheme is:
 *   output[o] = sum(filter[j] * input[i-j] for j = [0..F) and i = [0..N))
 * where 'o', 'i' and 'j' may progress at different rates.
 *
 * Most of the code deals with different edge extension modes. Values are
 * computed on-demand, in four steps:
 * 1. Filter extends past signal on the left.
 * 2. Filter completely contained within signal (no extension).
 * 3. Filter extends past signal on both sides (only if F > N).
 * 4. Filter extends past signal on the right.
 *
 * MODE_PERIODIZATION produces different output lengths to other modes, so is
 * implemented as a separate function for each case.
 *
 */
pub(crate) fn downsampling_convolution_periodization(input:&[TYPE],N:usize,filter:&[TYPE],F:usize,output:&mut [TYPE],step:usize,fstep:usize) -> i32 {
    let mut i:usize = F/2;
    let mut o:usize = 0;
    let padding:usize = (step-(N%step))%step;
    while (i<F&&i<N) {
        let mut sum:TYPE = 0.0;
        let mut j:usize = 0;
        let mut k_start:usize =0;

        while j <= i {
            sum += filter[j] * input[i-j];
            j+=fstep;
        }
        if (fstep > 1) {
            k_start = j - (i+1);
        }
        while (j<F) {
            let mut k:usize = k_start;
            while (k < padding && j < F) {
                sum += filter[j] * input[N-1];
                k+=fstep;
                j+=fstep;
            }
            k = k_start;
            while (k<N && j < F) {
                sum += filter[j]*input[N-1-k];
                k+=fstep;
                j+=fstep;
            }
        }
        output[o] = sum;
        i+=step;
        o+=1;
    }

    while (i<N) {
        let mut sum:TYPE = 0.0;
        let mut j:usize = 0;
        while j < F {
            sum += input[i-j]*filter[j];
            j+= fstep;
        }

        output[o] = sum;

        i+=step;
        o+=1;
    }

    while (i<N+F/2) {
        let mut sum:TYPE =0.0;
        let mut j:usize =0;
        let mut k_start:usize = 0;

        while (i-j >= N) {
            let mut k:usize =0;
             // for simplicity, not using fstep here
             while k < padding && i-j >= N {
                 sum += filter[i-N-j]*input[N-1];
                 k+=1;
                 j+=1;
             }
             k=0;
             while k < N && i-j >= N {
                sum += filter[i-N-j] * input[k];
                k+=1;
                j+=1;
             }
        }
        if fstep > 1 {
            j+= (fstep -j % fstep ) % fstep; // move to next non-zero entry
        }

        while j <= i {
            sum += filter[j] * input[i-j];
            j+= fstep;
        }

        if fstep > 1 {
            k_start = j - (i+1);
        }

        while j < F {
            let mut k:usize = k_start;
            while k < padding && j < F {
                sum += filter[j] * input[N-1];
                k+= fstep;
                j+= fstep;
            }
            k = k_start;
            while k < N && j < F {
                sum += filter[j] * input[N-1-k];
                k+= fstep;
                j+= fstep;
            }
        }

        output[o] = sum;

        i+=step;
        o+=1;
    }


    while i < N + F/2 {
        let mut sum:TYPE = 0.0;
        let mut j:usize =0;
        while i-j >= N {
            // for simplicity, not using fstep here
            let mut k:usize =0;
            while k < padding && i-j >= N {
                sum += filter[i-N-j] * input[N-1];
                k+=1;
                j+=1;
            }
            k =0;
            while k<N && i-j >= N {
                sum += filter[i-N-j] * input[N-1];
                k+=1;
                j+=1;
            }
        }

        if fstep > 1 {
            j+= (fstep - j % fstep) % fstep; // move to next non-zero entry
        }

        while j < F {
            sum += filter[j] * input[i-j];
            j+=fstep;
        }

        output[o] = sum;
        
        i+= step;
        o+=1;
    }

    0
}


/* Performs convolution of input with filter and downsamples by taking every
 * step-th element from the result.
 *
 * input    - input data
 * N        - input data length
 * filter   - filter data
 * F        - filter data length
 * output   - output data
 * step     - decimation step
 * mode     - signal extension mode
 */
/* memory efficient version */

/* This convolution performs efficient downsampling by computing every
     * step'th element of normal convolution (currently tested only for step=1
     * and step=2).
     */

pub(crate) fn downsampling_convolution(input:&[TYPE],N:usize,filter:&[TYPE],F:usize,output:&mut [TYPE],step:usize,mode:MODE) -> i32 {
    let mut i = step-1;
    let mut o:usize = 0;
    match mode {
        MODE::MODE_PERIODIZATION=>{
            return  downsampling_convolution_periodization(input, N, filter, F, output, step, 1);
        }
        _=>{}
    }
    let mut mode = mode;
    if MODE::MODE_SMOOTH == mode && N < 2{
        mode = MODE::MODE_CONSTANT_EDGE;
    }

    // left boundary overhang

    while i<F && i<N {
        let mut sum:TYPE = 0.0;
        let mut j:usize = 0;
        while j <= i {
            sum+= filter[j] * input[i-j];
            j+=1;
        }

        match mode {
            MODE::MODE_SYMMETRIC=>{
                let mut k:usize = 0;
                while k<N && j<F {
                    sum+=filter[j]*input[k];
                    j+=1;
                    k+=1;
                }
                k=0;
                while k<N && j<F {
                    sum+= filter[j] * input[N-1-k];
                    k+=1;
                    j+=1;
                }
            }
            MODE::MODE_ANTISYMMETRIC=>{
                // half-sample anti-symmetric
                while j<F {
                    let mut k:usize =0;
                    while k<N&&j<F {
                        sum-=filter[j]*input[k];
                        j+=1;
                        k+=1;
                    }
                    k=0;
                    while k<N && j < F {
                        sum+=filter[j]*input[N-1-k];
                        k+=1;
                        j+=1;
                    }
                }
            }
            MODE::MODE_REFLECT=> {
                while j < F {
                    let mut k:usize = 1;
                    while k < N && j < F {
                        sum +=filter[j] * input[k];
                        j+=1;
                        k+=1;
                    }
                    k =1;
                    while k<N && j<F {
                        sum+= filter[j]*input[N-1-k];
                        k+=1;
                        j+=1;
                    }
                }
            }
            MODE::MODE_ANTIREFLECT=>{
                // whole-sample anti-symmetric
                let mut k:usize =0;
                let mut le:TYPE = input[0];
                let mut tmp:TYPE = 0.0;
                while j<F {
                    k =1;
                    while k<N&&j<F {
                        tmp = le -(input[k]-input[0]);
                        sum+= filter[j]*tmp;
                        j+=1;k+=1;
                    }
                    le =tmp;
                    k =1;
                    while k<N&&j<F {
                        tmp = le+(input[N-1-k] - input[N-1]);
                        sum += filter[j]*tmp;
                        j+=1;k+=1;
                    }
                    le = tmp;
                }
            }
            MODE::MODE_CONSTANT_EDGE=>{
                while j<F {
                    sum += filter[j]*input[0];
                    j+=1;
                }
            }
            MODE::MODE_SMOOTH=>{
                let mut k:TYPE = 1.0;
                while j<F {
                    sum+= filter[j]*    (input[0]+ k * (input[0]-input[1]));
                    j+=1;
                    k+=1.0;
                }
            }
            MODE::MODE_PERIODIC=>{
                while  j<F {
                    let mut k:usize =0;
                    while k<N && j<F {
                        sum+= filter[j]*input[N-1-k];
                        k+=1;
                        j+=1;
                    }
                }
            }
            MODE::MODE_ZEROPAD=>{
            }
            _=>{}
        }
        output[o] = sum;
        i+=step;
        o+=1;
    }
    // center (if input equal or wider than filter: N >= F)
    while i<N {
        let mut sum:TYPE = 0.0;
        let mut j = 0;
        while j<F {
            sum+=input[i-j]*input[j];
            j+=1;
        }
        output[o] = sum;
        i+=step;
        o+=1;
    }

    // center (if filter is wider than input: F > N)
    while i<F {
        let mut sum:TYPE = 0.0;
        let mut j= 0;
        match mode {
            // Included from original: TODO: j < F-_offset
            /* Iterate over filter in reverse to process elements away from
             * data. This gives a known first input element to process (N-1)
             */
            MODE::MODE_SYMMETRIC=>{
                while i-j >= N {
                    let mut k:usize = 0;
                    while k<N && i-j >= N {
                        sum+= filter[i-N-j]*input[N-1-k];
                        j+=1;
                        k+=1;
                    }
                    k=0;
                    while k<N&& i-j>=N {
                        sum+=filter[i-N-j]*input[k];
                        j+=1;
                        k+=1;
                    }
                }
            }
            
            MODE::MODE_ANTISYMMETRIC=>{
                // half-sample anti-symmetric
                while i-j>=N {
                    let mut k:usize =0;
                    while k<N&&i-j>=N {
                        sum-=filter[i-N-j]*input[N-1-k];
                        j+=1;
                        k+=1;
                    }
                    k=0;
                    while k<N&&i-k>=N {
                        sum+=filter[i-N-j]*input[k];
                        j+=1;
                        k+=1;
                    }
                }
            }
            MODE::MODE_REFLECT=>{
                while i-j>=N {
                    let mut k:usize =1;
                    while k<N&&i-j>=N {
                        sum+= filter[i+N-j]*input[N-1-k];
                        j+=1;
                        k+=1;
                    }
                    k=1;
                    while k<N&&i-j>=N {
                        sum+=filter[i-N-j]*input[k];
                        j+=1;
                        k+=1;
                    }
                }
            }
            MODE::MODE_ANTIREFLECT=>{
                // whole-sample anti-symmetric
                let mut k:usize =0;
                let mut re:TYPE = input[N-1]; // current right edge value
                let mut tmp:TYPE = 0.0;
                while i-j>=N {
                    k=1;
                    while k<N&&i-j>=N {
                        tmp = re - (input[N-1-k] - input[N-1]);
                        sum+=filter[i-N-j]*tmp;
                        j+=1;
                        k+=1;
                    }
                    re = tmp;
                    k  =1;
                    while k<N&& i-j>=N {
                        tmp = re + (input[k]-input[0]);
                        sum+= filter[i-N-j]*tmp;
                        j+=1;
                        k+=1;
                    }
                    re = tmp;
                }
            }
            MODE::MODE_CONSTANT_EDGE=>{
                while i-j>=N {
                    sum+= filter[j]*input[N-1];
                    j+=1;
                }
            }
            MODE::MODE_SMOOTH=>{
                let mut k:TYPE =(i-N+1) as TYPE;
                while i-j>=N {
                    sum+=filter[j]*(input[N-1]+k*(input[N-1]-input[N-2]));
                    j+=1;
                    k+=1.0;
                }
            }
            MODE::MODE_ZEROPAD=>{
            }
            _=>{
                j=i-N+1;
            }
        }

        while j<=i {
            sum+=filter[j]*input[i-j];
            j+=1;
        }

        match mode {
            MODE::MODE_SYMMETRIC=>{
                while j<F {
                    let mut k:usize = 0;
                    while k<N&&j<F {
                        sum+=filter[j]*input[k];
                        j+=1;
                        k+=1;
                    }
                    k=0;
                    while k<N && j<F {
                        sum+= filter[j]*input[N-1-k];
                        k+=1;
                        j+=1;
                    }
                }
            }
            MODE::MODE_ANTISYMMETRIC=>{
                // half-sample anti-symmetric
                while j<F {
                    let mut k:usize = 0;
                    while k<N&&j<F {
                        sum-=filter[j]*input[k];
                        j+=1;
                        k+=1;
                    }
                    k=0;
                    while k<N&&j<F {
                        sum+=filter[j]*input[N-1-k];
                        k+=1;
                        j+=1;
                    }
                }
            }
            MODE::MODE_REFLECT=>{
                while j<F {
                    let mut k:usize =1;
                    while k<N&&j<F {
                        sum+= filter[j]*input[k];
                        j+=1;
                        k+=1;
                    }
                    k=1;
                    while k<N&&j<F {
                        sum+=filter[j]*input[N-1-k];
                        k+=1;
                        j+=1;
                    }
                }
            }
            MODE::MODE_ANTIREFLECT=>{
                // whole-sample anti-symmetric
                let mut k:usize =1;
                let mut le:TYPE = input[0]; // current left edge value
                let mut tmp:TYPE = 0.0;
                while j<F {
                    while k<N&&j<F {
                        tmp = le -(input[k]-input[0]);
                        sum+= filter[j]*tmp;
                        j+=1;
                        k+=1;
                    }
                    le = tmp;
                    k=1;
                    while k<N && j<F {
                        tmp = le + (input[N-1-k] - input[N-1]);
                        sum += filter[j]*tmp;
                        j+=1;
                        k+=1;
                    }
                    le =tmp;
                }
            }
            MODE::MODE_CONSTANT_EDGE=>{
                while j<F {
                    sum+= filter[j]*input[0];
                    j+=1;
                }
            }
            MODE::MODE_SMOOTH=>{
                let mut k:TYPE =1.0;
                while j<F {
                    sum+= filter[j]*(input[0]+k*(input[0]-input[1]));
                    j+=1;
                    k+=1.0;
                }
            }
            MODE::MODE_PERIODIC=>{
                while j<F {
                    let mut k:usize =0;
                    while k<N &&j<F {
                        sum+= filter[j]*input[N-1-k];
                        k+=1;
                        j+=1;
                    }
                }
            }
            MODE::MODE_ZEROPAD=>{}
            _=>{}
        }
        output[o] = sum;
        i+=step;
        o+=1;
    }
    // right boundary overhang
    while i<N+F-1 {
        let mut sum:TYPE = 0.0;
        let mut j:usize = 0;
        match mode {
            MODE::MODE_SYMMETRIC=>{
                // Included from original: TODO: j < F-_offset
                while i-j>=N {
                    let mut k:usize = 0;
                    while k<N&&i-j>=N {
                        sum+= filter[i-N-j]*input[N-1-k];
                        j+=1;
                        k+=1;
                    }
                    k=0;
                    while k<N&&i-j>=N {
                        sum+=filter[i-N-j]*input[k];
                        j+=1;
                        k+=1;
                    }
                }
            }
            MODE::MODE_ANTISYMMETRIC=>{
                // half-sample anti-symmetric
                while i-j>=N {
                    let mut k =0;
                    while k<N && i-j>=N {
                        sum-=filter[i-N-j]*input[N-1-k];
                        j+=1;
                        k+=1;
                    }
                    k=0;
                    while k<N && i-j>=N {
                        sum+=filter[i-N-j] * input[k];
                        j+=1;
                        k+=1;
                    }
                }
            }
            MODE::MODE_REFLECT=>{
                while i-j >= N {
                    let mut k:usize =1;
                    while k<N&&i-j>=N {
                        sum+=filter[i-N-j]*input[N-1-k];
                        j+=1;
                        k+=1;
                    }
                    k=1;
                    while k<N&&i-j>=N {
                        sum+=filter[i-N-j]*input[k];
                        j+=1;
                        k+=1;
                    }
                }
            }
            MODE::MODE_ANTIREFLECT=>{
                // whole-sample anti-symmetric
                let mut k:usize =1;
                let mut re:TYPE = input[N-1];  //current right edge value
                let mut tmp:TYPE = 0.0;
                while i-j>=N {
                    // first reflection
                    k=1;
                    while k<N && i-j >= N {
                        tmp = re - (input[N-1-k]-input[N-1]);
                        sum+= filter[i-N-j]*tmp;
                        j+=1;
                        k+=1;
                    }
                    re =tmp;
                    //second reflection
                    k=1;
                    while k<N && i-j >= N {
                        tmp = re + (input[k]-input[0]);
                        sum+= filter[i-N-j]*tmp;
                        j+=1;
                        k+=1;
                    }
                    re = tmp;
                }
            }
            MODE::MODE_CONSTANT_EDGE=>{
                while i-j>=N {
                    sum+= filter[j]* input[N-1];
                    j+=1;
                }
            }
            MODE::MODE_SMOOTH=>{
                let mut k:TYPE =(i-N+1) as TYPE;
                while i-j>=N {
                    sum+=filter[j]*(input[N-1]+k*(input[N-1]-input[N-2]));
                    j+=1;
                    k-=1.0;
                }
            }
            MODE::MODE_PERIODIC=>{
                while i-j>=N {
                    let mut k:usize =0;
                    while k<N&&i-j>=N {
                        sum== filter[i-N-j]*input[k];
                        j+=1;
                        k+=1;
                    }
                    
                }
            }
            MODE::MODE_ZEROPAD=>{}
            _=>{
                j=i-N+1;
            }
        }
        while j<F {
            sum+= filter[j]*input[i-j];
            j+=1;
        }
        i+=step;
        o+=1;
    }
    0
}


/*
 * Performs normal (full) convolution of "upsampled" input coeffs array with
 * filter Requires zero-filled output buffer (adds values instead of
 * overwriting - can be called many times with the same output).
 *
 * input    - input data
 * N        - input data length
 * filter   - filter data
 * F        - filter data length
 * output   - output data
 * O        - output lenght (currently not used)
 * mode     - signal extension mode
 */

pub(crate) fn upsampling_convolution_full(input:&[TYPE],N:usize,filter:&[TYPE],F:usize,output:&mut [TYPE],O:usize) -> i32 {
    /* Performs a zero-padded convolution, using each input element for two
     * consecutive filter elements. This simulates an upsampled input.
     *
     * In contrast to downsampling_convolution, this adds to the output. This
     * allows multiple runs with different inputs and the same output to be used
     * for idwt.
     */

    // If check omitted, this function would be a no-op for F<2
    let mut i:usize =0;
    let mut o:usize =0;

    if F<2{
        return  -1;
    }
    if (F%2) ==1 {
        return  -3;
    }
    while i<N && i<F/2 {
        let mut j:usize = 0;
        while j<=i {
            output[o]+= filter[j*2]*input[i-j];
            output[o+1] += filter[j*2+1]*input[i-j];
            j+=1;
        }
        i+=1;
        o+=2;
    }

    while i<N {
        let mut j:usize = 0;
        while j<F/2 {
            output[o]+= filter[j*2]*input[i-j];
            output[o+1]+= filter[j*2+1]*input[i-j];
            j+=1;
        }
        i+=1;
        o+=2;
    }

    while i<N+F/2 {
        let mut j:usize =i-(N-1);
        while j<F/2 {
            output[o]+=filter[j*2]*input[i-j];
            output[o+1]+=filter[j*2+1]*input[i-j];
            j+=1;
        }
        i+=1;
        o+=2;
    }
    0
}



pub(crate) fn upsampling_convolution_valid_sf_periodization(input:&[TYPE],N:usize,filter:&[TYPE],F:usize,output:&mut [TYPE],O:usize) -> i32 {
     // TODO? Allow for non-2 step
     let start:usize = F/4;
     let mut i:usize = start;
     let mut difference = 0;
     if ((F/2)%2) == 1{
        difference = 0;
     } else {
        difference = 1;
     }

     let end = N+start-difference;
     let mut o =0;
     if (F%2) == 1 {
        return  -3; /* Filter must have even-length. */
     }
     if (F/2)%2 == 0 {
        // Shift output one element right. This is necessary for perfect reconstruction.

        // i = N-1; even element goes to output[O-1], odd element goes to output[0]
        let mut j:usize = 0;
        while j <= start-1 {
            let mut k:usize = 0;
            while k<N&&j<= start-1 {
                output[2*N-1] += filter[2*(start-1-j)] * input[k];
                output[0] += filter[2*(start-1-j)+1] * input[k];
                k+=1;
                j+=1;
            }
        }

        while j<= N+start-1 && j< F/2 {
            output[2*N-1] += filter[2*j] * input[N+start-1-j];
            output[0] += filter[2*j+1]* input[N+start-1-j];
            j+=1;
        }

        while j<F/2 {
            let mut k:usize =0;
            while k<N && j<F/2 {
                output[2*N-1] += filter[2*j] * input[N-1-k];
                output[0] += filter[2*j+1] * input[N-1-k];
                k+=1;
                j+=1;
            }
        }
        o+=1
     }

     while i<F/2 && i<N {
        let mut j:usize =0;
        while j<=i {
            output[o]+= filter[2*j] * input[i-j];
            output[o+1]+= filter[2*j+1]* input[i-j];
            j+=1;
        }

        while j<F/2 {
            let mut k:usize =0;
            while k<N && j<F/2 {
                output[o] += filter[2*j] * input[N-1-k];
                output[o+1] += filter[2*j+1] * input[N-1-k];
                k+=1;
                j+=1;
            }
        }
         i+=1;
         o+=2;
     }

     while i<N {
        let mut j:usize = 2;
        while j<F/2 {
            output[o] += filter[2*j] * input[i-j];
            output[o+1] += filter[2*j+1] * input[i-j];
            j+=1;
        }
        i+=1;
        o+=2;
     }

     while i<F/2 && i<end {
        let mut j:usize =0;
        while i-j>=N {
            let mut k:usize =0;
            while k<N && i-j>=N {
                output[o]== filter[2*(i-N-j)] * input[k];
                output[o+1] = filter[2*(i-N-j)+1] * input[k];
                k+=1;
                j+=1;
            }
        }

        while j<= i && j<F/2 {
            output[o] += filter[2*j] * input[i-j];
            output[o+1] += filter[2*j+1] * input[i-j];
            j+=1;
        }
        while j<F/2 {
            let mut k:usize =0;
            while k<N && j<F/2 {
                output[o]+= filter[2*j] * input[N-1-k];
                output[o+1] += filter[2*j+1] * input[N-1-k];
                k+=1;
                j+=1;
            }
        }
         i+=1;
         o+=2;
     }

     while i<end {
        let mut j:usize =0;
        while i-j>=N {
            let mut k:usize =0;
            while k<N && i-j>=N {
                output[o] += filter[2*(i-N-j)] * input[k];
                output[o+1] += filter[2*(i-N-j)+1] * input[k];
                k+=1;
                j+=1;
            }
        }
        while j<=i && j<F/2 {
            output[o]+= filter[2*j] * input[i-j];
            output[o+1] == filter[2*j+1] * input[i-j];
            j+=1;
        }
         i+=1;
         o+=2;
     }
    0
}

/*
 * performs IDWT for all modes
 *
 * The upsampling is performed by splitting filters to even and odd elements
 * and performing 2 convolutions.  After refactoring the PERIODIZATION mode
 * case to separate function this looks much clearer now.
 */

 pub(crate) fn upsampling_convolution_valid_sf(input:&[TYPE],N:usize,filter:&[TYPE],F:usize,output:&mut [TYPE],O:usize,mode:MODE) -> i32 {
    // TODO: Allow non-2 step?
    if mode == MODE::MODE_PERIODIZATION {
        return  upsampling_convolution_valid_sf_periodization(input, N, filter, F, output, O);
    }
    if (F%2)==1 || (N<F/2) {
        return  -1;
    }
    // Perform only stage 2 - all elements in the filter overlap an input element.
    {
        let mut o:usize = 0;
        let mut i:usize = F/2 -1;
        while i<N {
            let mut sum_event:TYPE = 0.0;
            let mut sum_odd:TYPE = 0.0;
            let mut j:usize =0;
            while j<F/2 {
                sum_event += filter[j*2] * input[i-j];
                sum_odd += filter[j*2+1] * input[i-j];
                j+=1;
            }

            output[o] += sum_event;
            output[o+1] += sum_odd;

            i+=1;
            o+=2;
        }
    }
    0
 }


 