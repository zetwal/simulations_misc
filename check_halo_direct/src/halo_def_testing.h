#include <string>

// List all the halo float fields here
//#define N_HALO_FLOATS 73
#define N_HALO_FLOATS 37
#define N_HALO_FLOATS_SOD 5 
#define N_HALO_INTS_SOD 1
#define N_HALO_FLOATS_E 18

// these are referenced in tree building -- all others are carried along
// must retain this order for these codes - do not change this line
enum named_fields {fof_halo_mass=0, fof_halo_center_x=2, fof_halo_center_y=3, fof_halo_center_z=4};

const std::string float_var_names_test[N_HALO_FLOATS] = {
"fof_halo_mass",
"fof_halo_ke",
"fof_halo_center_x",
"fof_halo_center_y",
"fof_halo_center_z",
"fof_halo_angmom_x",
"fof_halo_angmom_y",
"fof_halo_angmom_z",
"fof_halo_max_cir_vel",
"fof_halo_com_x",
"fof_halo_com_y",
"fof_halo_com_z",
"fof_halo_com_vx",
"fof_halo_com_vy",
"fof_halo_com_vz",
"fof_halo_1D_vel_disp",
"sod_halo_radius",
"sod_halo_mass",
"sod_halo_ke",
"sod_halo_1D_vel_disp",
"sod_halo_max_cir_vel",
"sod_halo_center_x",
"sod_halo_center_y",
"sod_halo_center_z",
"sod_halo_angmom_x",
"sod_halo_angmom_y",
"sod_halo_angmom_z",
"sod_halo_com_x",
"sod_halo_com_y",
"sod_halo_com_z",
"sod_halo_com_vx",
"sod_halo_com_vy",
"sod_halo_com_vz",
"sod_halo_cdelta",
"sod_halo_cdelta_error",
"sod_halo_c_acc_mass",
"sod_halo_c_peak_mass"
};





const std::string float_var_names_test2[N_HALO_FLOATS] = {
		"fof_halo_mass",
"fof_halo_ke",
"fof_halo_center_x",
"fof_halo_center_y",
"fof_halo_center_z",
"fof_halo_angmom_x",
"fof_halo_angmom_y",
"fof_halo_angmom_z",
"fof_halo_max_cir_vel",
"fof_halo_com_x",
"fof_halo_com_y",
"fof_halo_com_z",
"fof_halo_com_vx",
"fof_halo_com_vy",
"fof_halo_com_vz",
"fof_halo_1D_vel_disp",
"sod_halo_radius",
"sod_halo_mass",
"sod_halo_ke",
"sod_halo_1D_vel_disp",
"sod_halo_max_cir_vel",
"sod_halo_center_x",
"sod_halo_center_y",
"sod_halo_center_z",
"sod_halo_angmom_x",
"sod_halo_angmom_y",
"sod_halo_angmom_z",
"sod_halo_com_x",
"sod_halo_com_y",
"sod_halo_com_z",
"sod_halo_com_vx",
"sod_halo_com_vy",
"sod_halo_com_vz",
"sod_halo_cdelta",
"sod_halo_cdelta_error",
"sod_halo_c_acc_mass",
"sod_halo_c_peak_mass"
};



const std::string float_var_names_ellipticity[N_HALO_FLOATS_E] = {
        "sod_halo_eigS1X",
        "sod_halo_eigS1Y",
        "sod_halo_eigS1Z",
        "sod_halo_eigS2X",
        "sod_halo_eigS2Y",
        "sod_halo_eigS2Z",
        "sod_halo_eigS3X",
        "sod_halo_eigS3Y",
        "sod_halo_eigS3Z",
        "sod_halo_eigR1X",
        "sod_halo_eigR1Y",
        "sod_halo_eigR1Z",
        "sod_halo_eigR2X",
        "sod_halo_eigR2Y",
        "sod_halo_eigR2Z",
        "sod_halo_eigR3X",
        "sod_halo_eigR3Y",
        "sod_halo_eigR3Z"
};


const std::string float_var_names_sodbin[N_HALO_FLOATS_SOD] = {
	"sod_halo_bin_mass",
	"sod_halo_bin_radius",
	"sod_halo_bin_rho",
	"sod_halo_bin_rho_ratio",
	"sod_halo_bin_rad_vel"
};
const std::string int_var_names_sodbin[N_HALO_INTS_SOD] = {
	"sod_halo_bin_count"
};

