"""This module contains CUDA source code for the evaluation of the electric 
and magnetic field from a periodic VWF expansion."""


# This cuda kernel is used for the evaluation of the electric field of
# periodic plane wave expansions.
pwe_periodic_electric_field_evaluation_code = """
    #define LEN_X %i
    #define LEN_K %i
    #define LEN_A %i
    #define RE_ONE_OVER_K %.10f
    #define IM_ONE_OVER_K %.10f
    
    __device__ void inner_sum(const int i_k, const float x, const float y, const float z, 
                                   const float re_kp, const float im_kp, const float re_kz, const float im_kz, 
                                   const float *alpha_array, 
                                   const float *re_g_te_array, const float *im_g_te_array,
                                   const float *re_g_tm_array, const float *im_g_tm_array,
                                   float *re_summand_x, float *im_summand_x,
                                   float *re_summand_y, float *im_summand_y,
                                   float *re_summand_z, float *im_summand_z)
    {
        float a = 0.0;
        
        float re_x_result = 0.0;
        float im_x_result = 0.0;
        float re_y_result = 0.0;
        float im_y_result = 0.0;
        float re_z_result = 0.0;
        float im_z_result = 0.0;
        
        float re_e_x_inner_summand = 0.0;
        float im_e_x_inner_summand = 0.0;
        float re_e_y_inner_summand = 0.0;
        float im_e_y_inner_summand = 0.0;
        float re_e_z_inner_summand = 0.0;
        float im_e_z_inner_summand = 0.0;

        for (int i_a=0; i_a<LEN_A; i_a++)
        {
            a = alpha_array[i_a];
            float sin_a = sinf(a);
            float cos_a = cosf(a);
            
            float re_ikr = -im_kp * (cos_a * x + sin_a * y) - im_kz * z;
            float im_ikr = re_kp * (cos_a * x + sin_a * y) + re_kz * z;
            
            float ereikr = expf(re_ikr);
            float re_eikr =  ereikr * cosf(im_ikr);
            float im_eikr =  ereikr * sinf(im_ikr);
            
            int i_ka = i_k * LEN_A + i_a;
            
            // pol=TE
            float re_g = re_g_te_array[i_ka];
            float im_g = im_g_te_array[i_ka];

            float re_geikr = re_g * re_eikr - im_g * im_eikr;
            float im_geikr = re_g * im_eikr + im_g * re_eikr;
            
            re_e_x_inner_summand = -sin_a * re_geikr;
            im_e_x_inner_summand = -sin_a * im_geikr;
            re_e_y_inner_summand = cos_a * re_geikr;
            im_e_y_inner_summand = cos_a * im_geikr;

            // pol=TM
            re_g = re_g_tm_array[i_ka];
            im_g = im_g_tm_array[i_ka];

            re_geikr = re_g * re_eikr - im_g * im_eikr;
            im_geikr = re_g * im_eikr + im_g * re_eikr;
            
            float re_kzgeikr = re_kz * re_geikr - im_kz * im_geikr;
            float im_kzgeikr = re_kz * im_geikr + im_kz * re_geikr;
            
            float re_kzkgeikr = re_kzgeikr * RE_ONE_OVER_K - im_kzgeikr * IM_ONE_OVER_K;
            float im_kzkgeikr = re_kzgeikr * IM_ONE_OVER_K + im_kzgeikr * RE_ONE_OVER_K;
            
            re_e_x_inner_summand += cos_a * re_kzkgeikr;
            im_e_x_inner_summand += cos_a * im_kzkgeikr;
            re_e_y_inner_summand += sin_a * re_kzkgeikr;
            im_e_y_inner_summand += sin_a * im_kzkgeikr;

            float re_kpgeikr = re_kp * re_geikr - im_kp * im_geikr;
            float im_kpgeikr = re_kp * im_geikr + im_kp * re_geikr;
            re_e_z_inner_summand = -(re_kpgeikr * RE_ONE_OVER_K - im_kpgeikr * IM_ONE_OVER_K);
            im_e_z_inner_summand = -(re_kpgeikr * IM_ONE_OVER_K + im_kpgeikr * RE_ONE_OVER_K);
            

            re_x_result += re_e_x_inner_summand;
            im_x_result += im_e_x_inner_summand;
            re_y_result += re_e_y_inner_summand;
            im_y_result += im_e_y_inner_summand;
            re_z_result += re_e_z_inner_summand;
            im_z_result += im_e_z_inner_summand;

        }
        re_summand_x[0] = re_x_result;
        im_summand_x[0] = im_x_result;
        re_summand_y[0] = re_y_result;
        im_summand_y[0] = im_y_result;
        re_summand_z[0] = re_z_result;
        im_summand_z[0] = im_z_result;
    }
    
    
    __global__ void electric_field(const float *re_kp_array, const float *im_kp_array, 
                                   const float *re_kz_array, const float *im_kz_array,
                                   const float *alpha_array,  
                                   const float *x_array, const float *y_array, const float *z_array,
                                   const float *re_g_te_array, const float *im_g_te_array,
                                   const float *re_g_tm_array, const float *im_g_tm_array,
                                   float *re_e_x, float *im_e_x, float *re_e_y, float *im_e_y, 
                                   float *re_e_z, float *im_e_z)
    {
        unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= LEN_X) return;
        
        float x = x_array[i];
        float y = y_array[i];
        float z = z_array[i];
        
        float re_kp = 0.0;
        float im_kp = 0.0;

        float re_x_result = 0.0;
        float im_x_result = 0.0;
        float re_y_result = 0.0;
        float im_y_result = 0.0;
        float re_z_result = 0.0;
        float im_z_result = 0.0;
        
        for (int i_k=0; i_k<LEN_K; i_k++)
        {            
            re_kp = re_kp_array[i_k];
            im_kp = im_kp_array[i_k];

            float re_kz = re_kz_array[i_k];
            float im_kz = im_kz_array[i_k];
            
            float re_x_inner_sum = 0.0;
            float im_x_inner_sum = 0.0;
            float re_y_inner_sum = 0.0;
            float im_y_inner_sum = 0.0;
            float re_z_inner_sum = 0.0;
            float im_z_inner_sum = 0.0;
            
            inner_sum(i_k, x, y, z, re_kp, im_kp, re_kz, im_kz, alpha_array, 
                           re_g_te_array, im_g_te_array, re_g_tm_array, im_g_tm_array,
                           &re_x_inner_sum, &im_x_inner_sum, &re_y_inner_sum, &im_y_inner_sum, 
                           &re_z_inner_sum, &im_z_inner_sum);


            re_x_result += re_x_inner_sum;
            im_x_result += im_x_inner_sum;
            re_y_result += re_y_inner_sum;
            im_y_result += im_y_inner_sum;
            re_z_result += re_z_inner_sum;
            im_z_result += im_z_inner_sum;

        }
        re_e_x[i] = re_x_result;
        im_e_x[i] = im_x_result;
        re_e_y[i] = re_y_result;
        im_e_y[i] = im_y_result;
        re_e_z[i] = re_z_result;
        im_e_z[i] = im_z_result;
    }
"""




# This cuda kernel is used for the evaluation of the magnetic field of
# periodic plane wave expansions.
pwe_periodic_magnetic_field_evaluation_code = """
    #define LEN_X %i
    #define LEN_K %i
    #define LEN_A %i
    #define RE_K %.10f
    #define IM_K %.10f
    
    __device__ void inner_sum(const int i_k, const float x, const float y, const float z, 
                                   const float re_kp, const float im_kp, const float re_kz, const float im_kz, 
                                   const float *alpha_array, 
                                   const float *re_g_te_array, const float *im_g_te_array,
                                   const float *re_g_tm_array, const float *im_g_tm_array,
                                   float *re_summand_x, float *im_summand_x,
                                   float *re_summand_y, float *im_summand_y,
                                   float *re_summand_z, float *im_summand_z)
    {
        float a = 0.0;
        
        float re_x_result = 0.0;
        float im_x_result = 0.0;
        float re_y_result = 0.0;
        float im_y_result = 0.0;
        float re_z_result = 0.0;
        float im_z_result = 0.0;
        
        float re_h_x_inner_summand = 0.0;
        float im_h_x_inner_summand = 0.0;
        float re_h_y_inner_summand = 0.0;
        float im_h_y_inner_summand = 0.0;
        float re_h_z_inner_summand = 0.0;
        float im_h_z_inner_summand = 0.0;

        for (int i_a=0; i_a<LEN_A; i_a++)
        {
            a = alpha_array[i_a];
            float sin_a = sinf(a);
            float cos_a = cosf(a);
            
            float re_ikr = -im_kp * (cos_a * x + sin_a * y) - im_kz * z;
            float im_ikr = re_kp * (cos_a * x + sin_a * y) + re_kz * z;
            
            float ereikr = expf(re_ikr);
            float re_eikr =  ereikr * cosf(im_ikr);
            float im_eikr =  ereikr * sinf(im_ikr);
            
            int i_ka = i_k * LEN_A + i_a;
            
            // pol=TE
            float re_g = re_g_te_array[i_ka];
            float im_g = im_g_te_array[i_ka];

            float re_geikr = re_g * re_eikr - im_g * im_eikr;
            float im_geikr = re_g * im_eikr + im_g * re_eikr;
            
            float re_kzgeikr = re_kz * re_geikr - im_kz * im_geikr;
            float im_kzgeikr = re_kz * im_geikr + im_kz * re_geikr;
            
            re_h_x_inner_summand = -cos_a * re_kzgeikr;
            im_h_x_inner_summand = -cos_a * im_kzgeikr;
            re_h_y_inner_summand = -sin_a * re_kzgeikr;
            im_h_y_inner_summand = -sin_a * im_kzgeikr;
            
            re_h_z_inner_summand = re_kp * re_geikr - im_kp * im_geikr;
            im_h_z_inner_summand = re_kp * im_geikr + im_kp * re_geikr;
            
            // pol=TM
            re_g = re_g_tm_array[i_ka];
            im_g = im_g_tm_array[i_ka];

            re_geikr = re_g * re_eikr - im_g * im_eikr;
            im_geikr = re_g * im_eikr + im_g * re_eikr;
            
            float re_kgeikr = RE_K * re_geikr - IM_K * im_geikr;
            float im_kgeikr = RE_K * im_geikr + IM_K * re_geikr;
            
            re_h_x_inner_summand += -sin_a * re_kgeikr;
            im_h_x_inner_summand += -sin_a * im_kgeikr;
            re_h_y_inner_summand += cos_a * re_kgeikr;
            im_h_y_inner_summand += cos_a * im_kgeikr;
            
            re_x_result += re_h_x_inner_summand;
            im_x_result += im_h_x_inner_summand;
            re_y_result += re_h_y_inner_summand;
            im_y_result += im_h_y_inner_summand;
            re_z_result += re_h_z_inner_summand;
            im_z_result += im_h_z_inner_summand;

        }
        re_summand_x[0] = re_x_result;
        im_summand_x[0] = im_x_result;
        re_summand_y[0] = re_y_result;
        im_summand_y[0] = im_y_result;
        re_summand_z[0] = re_z_result;
        im_summand_z[0] = im_z_result;
    }
    
    
    __global__ void magnetic_field(const float *re_kp_array, const float *im_kp_array, 
                                   const float *re_kz_array, const float *im_kz_array,
                                   const float *alpha_array,  
                                   const float *x_array, const float *y_array, const float *z_array,
                                   const float *re_g_te_array, const float *im_g_te_array,
                                   const float *re_g_tm_array, const float *im_g_tm_array,
                                   float *re_h_x, float *im_h_x, float *re_h_y, float *im_h_y, 
                                   float *re_h_z, float *im_h_z)
    {
        unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= LEN_X) return;
        
        float x = x_array[i];
        float y = y_array[i];
        float z = z_array[i];
        
        float re_kp = 0.0;
        float im_kp = 0.0;

        float re_x_result = 0.0;
        float im_x_result = 0.0;
        float re_y_result = 0.0;
        float im_y_result = 0.0;
        float re_z_result = 0.0;
        float im_z_result = 0.0;
        
        for (int i_k=0; i_k<LEN_K; i_k++)
        {            
            re_kp = re_kp_array[i_k];
            im_kp = im_kp_array[i_k];

            float re_kz = re_kz_array[i_k];
            float im_kz = im_kz_array[i_k];
            
            float re_x_inner_sum = 0.0;
            float im_x_inner_sum = 0.0;
            float re_y_inner_sum = 0.0;
            float im_y_inner_sum = 0.0;
            float re_z_inner_sum = 0.0;
            float im_z_inner_sum = 0.0;
            
            inner_sum(i_k, x, y, z, re_kp, im_kp, re_kz, im_kz, alpha_array, 
                           re_g_te_array, im_g_te_array, re_g_tm_array, im_g_tm_array,
                           &re_x_inner_sum, &im_x_inner_sum, &re_y_inner_sum, &im_y_inner_sum, 
                           &re_z_inner_sum, &im_z_inner_sum);

            re_x_result += re_x_inner_sum;
            im_x_result += im_x_inner_sum;
            re_y_result += re_y_inner_sum;
            im_y_result += im_y_inner_sum;
            re_z_result += re_z_inner_sum;
            im_z_result += im_z_inner_sum;

        }
        re_h_x[i] = re_x_result;
        im_h_x[i] = im_x_result;
        re_h_y[i] = re_y_result;
        im_h_y[i] = im_y_result;
        re_h_z[i] = re_z_result;
        im_h_z[i] = im_z_result;
    }
"""
