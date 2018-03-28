typedef struct node_type {
  char  *name;
  int   x, y;
  int   num_connectnions;
} NODE;

NODE nodes[] = {
   {"/webmod/basin_topg.f", 59, 181, 0},
   {"/webmod/io_chem.f", 57, 248, 0},
   {"/webmod/obs_chem.f", 57, 372, 0},
   {"/webmod/obs_webmod.f", 57, 314, 0},
   {"/webmod/soltab_prms.f", 88, 16, 0},
   {"/webmod/ccsolrad_web.f", 240, 17, 4},
   {"/webmod/temp_1sta_prms.f", 154, 98, 1},
   {"/webmod/precip_web.f", 156, 166, 2},
   {"/webmod/irrig_web.f", 318, 166, 5},
   {"/webmod/intcp_prms.f", 242, 245, 7},
   {"/webmod/nwsmelt_topg.f", 244, 320, 2},
   {"/webmod/potet_hamon_prms.f", 319, 99, 2},
   {"/webmod/topmod_chem.f", 245, 392, 7},
   {"/webmod/top2clark.f", 248,463, 1},
   {"/webmod/route_clark.f", 249, 543, 4},
   {"/webmod/webmod_res.f", 436, 204, 9},
   {"/webmod/phreeq_mms.f", 437, 305, 9},
   {"/webmod/web_sum.f", 438, 408, 5},
   NULL
};

int  con_index[] = {
   3, 4, 6, 7,
   3,
   3, 6,
   3, 6, 7, 12, 14,
   0, 3, 5, 6, 7, 10, 11,
   6, 9,
   5, 6,
   1, 3, 7, 8, 9, 10, 11,
   12,
   1, 3, 8, 13,
   1, 7, 8, 9, 10, 11, 12, 13, 14,
   1, 2, 7, 8, 10, 11, 12, 14, 15,
   1, 5, 6, 11, 15
};
