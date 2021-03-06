#define m3dc1_scorec_init m3dc1_domain_init_
#define m3dc1_scorec_finalize m3dc1_domain_finalize_
#define m3dc1_plane_setnum m3dc1_plane_setnum_
#define m3dc1_plane_getnum m3dc1_plane_getnum_
#define m3dc1_plane_getid m3dc1_plane_getid_
#define m3dc1_plane_setphirange m3dc1_plane_setphirange_
#define m3dc1_plane_setphi m3dc1_plane_setphi_
#define m3dc1_plane_getphi m3dc1_plane_getphi_
#define m3dc1_model_load m3dc1_model_load_
#define m3dc1_model_print m3dc1_model_print_
#define m3dc1_model_setnumplane m3dc1_model_setnumplane_
#define m3dc1_model_getnumplane m3dc1_model_getnumplane_
#define m3dc1_model_getedge m3dc1_model_getedge_
#define m3dc1_model_setpbc m3dc1_model_setpbc_
#define m3dc1_model_getpbc m3dc1_model_getpbc_
#define m3dc1_model_getmincoord m3dc1_model_getmincoord_
#define m3dc1_model_getmaxcoord m3dc1_model_getmaxcoord_
#define m3dc1_mesh_load m3dc1_mesh_load_
#define m3dc1_mesh_build3d m3dc1_mesh_build3d_
#define m3dc1_mesh_setcoord m3dc1_mesh_setcoord_
#define m3dc1_mesh_getcoord m3dc1_mesh_getcoord_
#define m3dc1_mesh_getnument m3dc1_mesh_getnument_
#define m3dc1_mesh_getnumownent m3dc1_mesh_getnumownent_
#define m3dc1_mesh_getnumglobalent m3dc1_mesh_getnumglobalent_
#define m3dc1_mesh_setordering m3dc1_mesh_setordering_
#define m3dc1_mesh_getordering m3dc1_mesh_getordering_
#define set_adapt_smooth_factor set_adapt_smooth_factor_
#define m3dc1_mesh_adapt m3dc1_mesh_adapt_
#define m3dc1_node_getglobalid m3dc1_node_getglobalid_
#define m3dc1_ent_getgeomclass m3dc1_ent_getgeomclass_
#define m3dc1_ent_getadj m3dc1_ent_getadj_
#define m3dc1_ent_getnumadj m3dc1_ent_getnumadj_
#define m3dc1_ent_getownpartid m3dc1_ent_getownpartid_
#define m3dc1_ent_ismine m3dc1_ent_ismine_
#define m3dc1_node_getcoord m3dc1_node_getcoord_
#define m3dc1_node_getnormvec m3dc1_node_getnormvec_
#define m3dc1_node_isongeombdry m3dc1_node_isongeombdry_
#define m3dc1_field_create m3dc1_field_create_
#define m3dc1_field_delete m3dc1_field_delete_
#define m3dc1_field_exist m3dc1_field_exist_
#define m3dc1_field_sync m3dc1_field_sync_
#define m3dc1_field_sum m3dc1_field_sum_
#define m3dc1_ent_getlocaldofid m3dc1_ent_getlocaldofid_
#define m3dc1_ent_getglobaldofid m3dc1_ent_getglobaldofid_
#define m3dc1_ent_getnumdof m3dc1_ent_getnumdof_
#define m3dc1_ent_setdofdata m3dc1_ent_setdofdata_
#define m3dc1_ent_getdofdata m3dc1_ent_getdofdata_
#define m3dc1_field_getlocaldofid m3dc1_field_getlocaldofid_
#define m3dc1_field_getowndofid m3dc1_field_getowndofid_
#define m3dc1_field_getglobaldofid m3dc1_field_getglobaldofid_
#define m3dc1_field_getnumlocaldof m3dc1_field_getnumlocaldof_
#define m3dc1_field_getnumowndof m3dc1_field_getnumowndof_
#define m3dc1_field_getnumglobaldof m3dc1_field_getnumglobaldof_
#define m3dc1_field_getdataptr m3dc1_field_getdataptr_
#define m3dc1_matrix_create m3dc1_matrix_create_
#define m3dc1_matrix_freeze m3dc1_matrix_freeze_
#define m3dc1_matrix_delete m3dc1_matrix_delete_
#define m3dc1_matrix_insert m3dc1_matrix_insert_
#define m3dc1_matrix_add m3dc1_matrix_add_
#define m3dc1_matrix_setbc m3dc1_matrix_setbc_
#define m3dc1_matrix_solve m3dc1_matrix_solve_
#define m3dc1_matrix_multiply m3dc1_matrix_multiply_
#define m3dc1_matrix_write m3dc1_matrix_write_
#define m3dc1_matrix_print m3dc1_matrix_print_
#define m3dc1_matrix_flush m3dc1_matrix_flush_
#define m3dc1_model_getplaneid m3dc1_model_getplaneid_
#define m3dc1_node_getcurv m3dc1_node_getcurv_
#define m3dc1_field_getnewid m3dc1_field_genid_
#define m3dc1_field_add m3dc1_field_add_
#define m3dc1_field_mult m3dc1_field_mult_
#define m3dc1_field_assign m3dc1_field_assign_
#define m3dc1_field_copy m3dc1_field_copy_
#define m3dc1_field_insert m3dc1_field_insert_
#define m3dc1_field_isnan m3dc1_field_isnan_
#define m3dc1_matrix_getiternum m3dc1_matrix_getiternum_
#define m3dc1_matrix_insertblock m3dc1_matrix_insertblock_
#define m3dc1_field_sumsq m3dc1_field_sumsq_
#define m3dc1_field_compare m3dc1_field_compare_
#define m3dc1_field_write m3dc1_field_write_
#define m3dc1_field_print m3dc1_field_print_
#define m3dc1_matrix_setlaplacebc m3dc1_matrix_setlaplacebc_
#define m3dc1_field_retrieve m3dc1_field_retrieve_
#define m3dc1_field_set m3dc1_field_set_
#define m3dc1_region_getoriginalface m3dc1_region_getoriginalface_
#define m3dc1_matrix_setassembleoption m3dc1_matrix_setassembleoption_
#define m3dc1_field_sum_plane m3dc1_field_sum_plane_
#define m3dc1_field_write_local m3dc1_field_write_local_
#define adapt_by_field adapt_by_field_
#define m3dc1_field_printcompnorm m3dc1_field_printcompnorm_
#define m3dc1_mesh_write m3dc1_mesh_write_
#define output_face_data output_face_data_
#define adapt_by_error_field adapt_by_error_field_
#define sum_edge_data sum_edge_data_
#define get_node_error_from_elm get_node_error_from_elm_
#define m3dc1_field_max m3dc1_field_max_
#define set_mesh_size_bound set_mesh_size_bound_
#define set_adapt_p set_adapt_p_
#define m3dc1_epetra_create m3dc1_epetra_create_
#define m3dc1_epetra_delete m3dc1_epetra_delete_
#define m3dc1_epetra_insert m3dc1_epetra_insert_
#define m3dc1_epetra_addblock m3dc1_epetra_addblock_
#define m3dc1_epetra_setbc m3dc1_epetra_setbc_
#define m3dc1_epetra_setlaplacebc m3dc1_epetra_setlaplacebc_
#define m3dc1_epetra_multiply m3dc1_epetra_multiply_
#define m3dc1_epetra_freeze m3dc1_epetra_freeze_
#define m3dc1_epetra_print m3dc1_epetra_print_
#define m3dc1_solver_aztec m3dc1_solver_aztec_
#define m3dc1_solver_amesos m3dc1_solver_amesos_
#define m3dc1_solver_getnumiter m3dc1_solver_getnumiter_
#define m3dc1_ghost_load m3dc1_ghost_load_
#define m3dc1_ghost_delete m3dc1_ghost_delete_




