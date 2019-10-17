#ifndef CNV_H 
#define CNV_H


const int64_t FIX_LENGTH = (int64_t)(1e10);
const int LOCAL_LINE_MAX = 1<<20;
const int MAX_DEPTH      = 10000 * 100;

const int HAP_UNK     = 1 << 0;
const int HAP1_DEL    = 1 << 1;
const int HAP2_DEL    = 1 << 2;
const int HAP_ALL_DEL = 1 << 3;

const int HAP1_DUP    = 1 << (1+16);
const int HAP2_DUP    = 1 << (2+16);
const int HAP_ALL_DUP = 1 << (3+16);
const int flag_dup    = HAP1_DUP | HAP2_DUP | HAP_ALL_DUP;
const int flag_del    = HAP1_DEL | HAP2_DEL | HAP_ALL_DEL;
const int MAX_FRM_DEPTH  = 10000;
const int min_pass_score = 20;

const int error_level = 2; // 0 for none, 1 for normal, 2 for debug

class MergedCNVInterval
{
public:
    int32_t tid;
    int32_t start_idx;
    int32_t end_idx;
	int32_t prev_idx;
    int32_t next_idx;
    double mean_depth;
    double variance;
    int genotype;
    double gt_score;
    double distance_to_next;
    bool is_active;
    int8_t cnv_type;
    std::string filter;
    std::string sv_id;

    MergedCNVInterval()
    {
        sv_id = ".";
    }

    int32_t length(int bin_size) const 
    {
        return (end_idx - start_idx) * bin_size;
    }

    double calculate_mean_depth (const std::vector <std::vector <double>> & wg_2d_depth_vector)
    {
        mean_depth = 0.0;
       
        for (int32_t idx = start_idx; idx < end_idx; idx++)
        {
            mean_depth += wg_2d_depth_vector[tid][idx];
        }

        mean_depth /= (end_idx - start_idx);

        return mean_depth;
    }

    double calculate_variance (const std::vector <std::vector <double>> & wg_2d_depth_vector)
    {
        double error; 
        variance = 0.0;
        for (int32_t idx = start_idx; idx < end_idx; idx++)
        {
            error = wg_2d_depth_vector[tid][idx] - mean_depth;
            variance += error * error;
        }
 
        variance /= (end_idx - start_idx);
       
        return variance;
    }

    double calculate_sd () const
    {
        return sqrt(variance);
    }

    int calculate_genotype(const std::vector <std::vector <double>> & wg_2d_depth_vector, double wg_mean_depth, double wg_sd_mean, int chrom_copy_number)
    {
        double d[5];
        double prior[5] = {1e-4, 1e-2, 0.98, 1e-2, 1e-4};
        double llr[5];
        
        double llr_max;
        double z, z_norm;
        int n = end_idx - start_idx;

        if (chrom_copy_number == 1)
        {
            prior[0] = 1e-2;
            prior[1] = 0.98;
            prior[2] = 1e-2;
            prior[3] = 1e-4;
            prior[4] = 1e-5;

        }else if (chrom_copy_number == 0)
        {
            prior[0] = 0.98;
            prior[1] = 1e-2;
            prior[2] = 1e-4;
            prior[3] = 1e-5;
            prior[4] = 1e-6;

        }

        z_norm = (mean_depth - wg_mean_depth) * sqrt(n) / wg_sd_mean;

        for (int i = 0; i < 5; i++)
        {
            d[i] = i * wg_mean_depth / 2;
            z = (mean_depth - d[i]) * sqrt(n) / wg_sd_mean;
            llr[i] = (z_norm * z_norm - z * z)/ 2.0 + log(prior[i]) * log(n) * 2;
            llr[i] /= log(10);
        }

        llr_max = llr[0]; 
        genotype = 0;
        for (int i = 1; i < 5; i++)
        {
            if (llr[i] > llr_max)
            {
                llr_max = llr[i];
                genotype = i;
            }
        }
        gt_score = llr[genotype];
        cnv_type = genotype - chrom_copy_number;
        
        return genotype;
    }

    int output_bedpe(const CHR_INFO * chr_info, int bin_size, FILE * fp) const
    {
        char * ctg;
        int pos1, pos2; 
        std::string sv_type;

        ctg = chr_info->chrname_list->data_list[tid];
        pos1 = start_idx * bin_size;
        pos2 = end_idx * bin_size;
        
        if (cnv_type < 0)
        {
            sv_type = "DEL"; 
        }else if (cnv_type > 0){
            sv_type = "DUP";
        }else{
            sv_type = "CNV_NEUTRAL";
        }

        fprintf(fp, "%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t", ctg, pos1, pos1+1, ctg, pos2, pos2+1, sv_type.c_str(), sv_id.c_str(), filter.c_str());
        fprintf(fp, "%d\t%.f\t%s\tSVMETHOD=READ_DEPTH;DEPTH=%.2f;COPY_NUMBER=%d\n", length(bin_size), gt_score, filter.c_str(), mean_depth, genotype);
        return 0;
    }

};

class CNVInterval
{
public:
    int32_t tid;
    int32_t start_pos;
    int32_t end_pos;
    float depth[4] = {0.0}; // depth each hap_type (0, 1, 2, 4) 0 for unknown, 3 for sum of 0,1,2
    float high_qual_depth[4] = {0.0}; // depth of reads with high map qual
    int32_t frm_depth[4];

    CNVInterval()
    {
        tid = -1;
        start_pos = -1;
        end_pos = -1;
        depth[0] = depth[1] = depth[2] = depth[3] = 0.0;
        high_qual_depth[0] = high_qual_depth[1] = high_qual_depth[2] = high_qual_depth[3] = 0.0;
        frm_depth[0] = frm_depth[1] = frm_depth[2] = frm_depth[3] = 0.0;
    };

    int length() const
    {
        return end_pos - start_pos;
    }

    Interval to_interval(const CHR_INFO * chr_info) const
    {
        Interval itv;
        itv.tid = tid;
        itv.start_pos = start_pos;
        itv.end_pos = end_pos;
        itv.ctg = chr_info->chrname_list->data_list[tid];
        return itv;
    }

    float depth_rate(int hap_type) const
    {
        if (frm_depth[hap_type] > 0)
        {
            return depth[hap_type] / frm_depth[hap_type];
        }else{
            return 0.0;
        }
    }
};


class CNVCall : public CNVInterval 
{
public:
    float  score[4] = {0}; // variant score of the CNV call. 0, 1, and 2 for hap_type 0, 1 and 2. 3 for all hap_type. 
    uint32_t flag = 0;
    int8_t cnv_hap_type; 
    int8_t cnv_type; // 1 for dup, -1 for del
    int16_t num_pe_supp[4]; // number of pe support for each hap_type (3 for total)
    int32_t num_cnv_bin[4];
    std::string aux_info;
    std::string sv_id;

    CNVCall()
    {
        tid = -1;
        start_pos = -1;
        end_pos = -1;
        depth[0] = depth[1] = depth[2] = depth[3] = 0.0;
        high_qual_depth[0] = high_qual_depth[1] = high_qual_depth[2] = high_qual_depth[3] = 0.0;
        frm_depth[0] = frm_depth[1] = frm_depth[2] = frm_depth[3] = 0.0;
        score[0] = score[1] = score[2] = score[3] = 0.0;
        flag = 0;
        cnv_hap_type = -1;
        cnv_type = 0;
        num_pe_supp[0] = num_pe_supp[1] = num_pe_supp[2] = num_pe_supp[3] = 0;
        num_cnv_bin[0] = num_cnv_bin[1] = num_cnv_bin[2] = num_cnv_bin[3] = 0;
        aux_info = "";
        sv_id = ".";
    }
    
    int output_line(FILE * out_fp, const CHR_INFO * const chr_info, int format = 1) const
    {
        const char * chrname = tid2chrname(tid, chr_info);
        
        std::string sv_type, flt;
        if (cnv_type == -1)
        {
            sv_type = "DEL";
        }else if (cnv_type == 1){
            sv_type = "DUP";
        }else{
            sv_type = "UNK";
        }

        flt = "LowQual";
        if (cnv_type != 0 && score[cnv_hap_type] >= min_pass_score)
        {
            flt = "PASS";
        }
        if (format == 1)
        {
            fprintf(out_fp, "%s\t%d\t%d\t", chrname, start_pos,  end_pos);
        }else{
            fprintf(out_fp, "%s\t%d\t%d\t%s\t%d\t%d\t", chrname, start_pos, start_pos+1, chrname, end_pos, end_pos+1);
        }
        
        fprintf(out_fp, "%s\t%s\t%d\t%.f\t%s\t", sv_type.c_str(), sv_id.c_str(), length(), score[cnv_hap_type], flt.c_str());

        if (cnv_hap_type > 0)
        {
            std::string hap_type;
            if (cnv_hap_type == 1){
                hap_type = "hap1";
            }else if (cnv_hap_type == 2){
                hap_type = "hap2";
            }else if (cnv_hap_type == 3){
                hap_type = "both";
            } else { hap_type = "unknown"; }

            fprintf(out_fp, "SVMETHOD=READ_PAIR;SV_HAP_TYPE=%s;NUM_READ_PAIR=%d;DEPTH=%.2f", hap_type.c_str(), num_pe_supp[cnv_hap_type], depth[cnv_hap_type]);
        }else{
            fprintf(out_fp, "SVMETHOD=READ_PAIR;SV_HAP_TYPE=unknown");
        }
        
        return 0;
    }

    int output_bed_line(FILE * out_fp, const CHR_INFO * const chr_info) const
    {
        output_line(out_fp, chr_info, 1);
        return 0;
    }

    int output_bedpe_line(FILE * out_fp, const CHR_INFO * const chr_info) const
    {
        output_line(out_fp, chr_info, 2);
        return 0;
    }
};

class CNVSettings
{
public:
    const CHR_INFO * chr_info;
    int bin_size;
    double wg_mean_depth; 
    double wg_sd_depth;
    double alpha; // smooth parameter, very important for sensitivity. 
    int32_t min_large_cnv_size;
    double min_large_cnv_llr;
    std::map <int32_t, int32_t> special_chrom_tid_map;
    std::vector<int32_t> chrom_copy_number_vector;
    std::vector <Interval> gap_region_interval_vector;

    std::vector <QuantileNumbers> depth_quantiles_vector;
    std::vector <QuantileNumbers> hiqh_qual_depth_quantiles_vector;

    std::vector <QuantileNumbers> depth_rate_quantiles_vector;

    std::vector <float> upper_depth_cutoff;
    std::vector <float> lower_depth_cutoff;
    std::vector <float> upper_depth_pvalue;
    std::vector <float> lower_depth_pvalue;

    std::vector <float> upper_depth_ratio_cutoff;
    std::vector <float> lower_depth_ratio_cutoff;
    std::vector <float> upper_depth_ratio_pvalue;
    std::vector <float> lower_depth_ratio_pvalue;  

    
    CNVSettings()
    {
        for (size_t i = 0; i < 4; i++)
        {
            upper_depth_cutoff.push_back(0.0);
            lower_depth_cutoff.push_back(0.0);
            upper_depth_pvalue.push_back(0.0);
            lower_depth_pvalue.push_back(0.0);

            depth_quantiles_vector.push_back(QuantileNumbers());
            hiqh_qual_depth_quantiles_vector.push_back(QuantileNumbers());

            upper_depth_ratio_cutoff.push_back(0.0);
            lower_depth_ratio_cutoff.push_back(0.0);
            upper_depth_ratio_pvalue.push_back(0.0);
            lower_depth_ratio_pvalue.push_back(0.0);

            depth_rate_quantiles_vector.push_back(QuantileNumbers());
        }
    }
};


int generate_whole_genome_interval_vector(std::vector <std::vector <CNVInterval>> & wg_2d_interval_vector, const CHR_INFO * chr_info, int bin_size);

int calculate_depth_distributions(const std::vector <std::vector <CNVInterval> > & wg_2d_interval_vector, std::vector <QuantileNumbers> & out_depth_quantiles, std::vector <QuantileNumbers> & out_high_qual_depth_quantiles);

int get_depth_for_each_interval(std::string hap_type_read_depth_file, std::vector <std::vector <CNVInterval> > & wg_2d_interval_vector, const CHR_INFO * chr_info, int bin_size);


#endif