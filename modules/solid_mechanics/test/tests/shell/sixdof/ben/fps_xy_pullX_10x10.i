[Mesh]
#   [mesh]
#     type = FileMeshGenerator
#     file = flatplates_xy.e
#   []
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 10
    ny = 10
    xmin = 0.0
    xmax = 1.0
    ymin = 0.0
    ymax = 1.0
    show_info = true
  []
  [all_nodes]
    type = BoundingBoxNodeSetGenerator
    input = gmg
    bottom_left = '-1e6 -1e6 -1e6'
    top_right = '1e6 1e6 1e6'
    new_boundary = 'all_nodes'
  []
[]

[Variables]
  [disp_x]
    order = FIRST
    family = LAGRANGE
  []
  [disp_y]
    order = FIRST
    family = LAGRANGE
  []
  [disp_z]
    order = FIRST
    family = LAGRANGE
  []
  [rot_x]
    order = FIRST
    family = LAGRANGE
  []
  [rot_y]
    order = FIRST
    family = LAGRANGE
  []
  [rot_z]
    order = FIRST
    family = LAGRANGE
  []
[]

[AuxVariables]
  [react_disp_x]
  []
  [react_disp_y]
  []
  [react_disp_z]
  []
  [react_rot_x]
  []
  [react_rot_y]
  []
  [react_rot_z]
  []
  [stress_xx]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_yy]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_zz]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_xy]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_yx]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_yz]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_zy]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_xz]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_zx]
    order = CONSTANT
    family = MONOMIAL
  []
  [strain_xx]
    order = CONSTANT
    family = MONOMIAL
  []
  [strain_yy]
    order = CONSTANT
    family = MONOMIAL
  []
  [strain_zz]
    order = CONSTANT
    family = MONOMIAL
  []
  [strain_xy]
    order = CONSTANT
    family = MONOMIAL
  []
  [strain_yx]
    order = CONSTANT
    family = MONOMIAL
  []
  [strain_xz]
    order = CONSTANT
    family = MONOMIAL
  []
  [strain_zx]
    order = CONSTANT
    family = MONOMIAL
  []
  [strain_yz]
    order = CONSTANT
    family = MONOMIAL
  []
  [strain_zy]
    order = CONSTANT
    family = MONOMIAL
  []
  
[]

[BCs]
  [fix_x]
    type = DirichletBC
    variable = disp_x
    boundary = '3'
    value = 0.0
  []
  [fix_y]
    type = DirichletBC
    variable = disp_y
    boundary = '0 2'
    # boundary = all_nodes
    value = 0.0
  []
  [fix_z]
    type = DirichletBC
    variable = disp_z
    boundary = '3'
    # boundary = all_nodes
    value = 0.0
  []
  [fix_rot_x]
    type = DirichletBC
    variable = rot_x
    # boundary = '0 1 2 3'
    boundary = all_nodes
    value = 0.0
  []
  [fix_rot_y]
    type = DirichletBC
    variable = rot_y
    # boundary = '3'
    boundary = all_nodes
    value = 0.0
  []
  [fix_rot_z]
    type = DirichletBC
    variable = rot_z
    #boundary = '0 1 2 3'
    boundary = all_nodes
    value = 0.0
  []
  [xy_pull_x]
    type = DirichletBC
    variable = disp_x
    boundary = '1'
    value = 0.01
  []
[]

#[DiracKernels]
#  [point1]
#    type = ConstantPointSource
#    variable = disp_x
#    point = '1 0 1'
#    value = -2.5 # P = 10
#  []
#[]

[NodalKernels]
#  [fx]
#    type = UserForcingFunctionNodalKernel
#    boundary = '1'
#    function = 10
#    variable = 'disp_x'
#  []
# [./constraint_z]
#   type = PenaltyDirichletNodalKernel
#   variable = rot_z
#   value = 0
#   penalty = 1e6
# []
[]

[AuxKernels]
  [stress_xx]
    type = ADRankTwoAux
    variable = stress_xx
    rank_two_tensor = global_stress_t_points_0
    index_i = 0
    index_j = 0
    execute_on = TIMESTEP_END
  []
  [strain_xx]
    type = ADRankTwoAux
    variable = strain_xx
    rank_two_tensor = t_points_0_total_strain
    index_i = 0
    index_j = 0
    execute_on = TIMESTEP_END
  []
  [stress_yy]
    type = ADRankTwoAux
    variable = stress_yy
    rank_two_tensor = global_stress_t_points_0
    index_i = 1
    index_j = 1
    execute_on = TIMESTEP_END
  []
  [strain_yy]
    type = ADRankTwoAux
    variable = strain_yy
    rank_two_tensor = t_points_0_total_strain
    index_i = 1
    index_j = 1
    execute_on = TIMESTEP_END
  []
  [stress_zz]
    type = ADRankTwoAux
    variable = stress_zz
    rank_two_tensor = global_stress_t_points_0
    index_i = 2
    index_j = 2
    execute_on = TIMESTEP_END
  []
  [strain_zz]
    type = ADRankTwoAux
    variable = strain_zz
    rank_two_tensor = t_points_0_total_strain
    index_i = 2
    index_j = 2
    execute_on = TIMESTEP_END
  []
  [stress_xy]
    type = ADRankTwoAux
    variable = stress_xy
    rank_two_tensor = global_stress_t_points_0
    index_i = 0
    index_j = 1
    execute_on = TIMESTEP_END
  []
  [strain_xy]
    type = ADRankTwoAux
    variable = strain_xy
    rank_two_tensor = t_points_0_total_strain
    index_i = 0
    index_j = 1
    execute_on = TIMESTEP_END
  []
  [stress_yx]
    type = ADRankTwoAux
    variable = stress_yx
    rank_two_tensor = global_stress_t_points_0
    index_i = 1
    index_j = 0
    execute_on = TIMESTEP_END
  []
  [strain_yx]
    type = ADRankTwoAux
    variable = strain_yx
    rank_two_tensor = t_points_0_total_strain
    index_i = 1
    index_j = 0
    execute_on = TIMESTEP_END
  []
  [stress_xz]
    type = ADRankTwoAux
    variable = stress_xz
    rank_two_tensor = global_stress_t_points_0
    index_i = 0
    index_j = 2
    execute_on = TIMESTEP_END
  []
  [strain_xz]
    type = ADRankTwoAux
    variable = strain_xz
    rank_two_tensor = t_points_0_total_strain
    index_i = 0
    index_j = 2
    execute_on = TIMESTEP_END
  []
  [stress_zx]
    type = ADRankTwoAux
    variable = stress_zx
    rank_two_tensor = global_stress_t_points_0
    index_i = 2
    index_j = 0
    execute_on = TIMESTEP_END
  []
  [strain_zx]
    type = ADRankTwoAux
    variable = strain_zx
    rank_two_tensor = t_points_0_total_strain
    index_i = 2
    index_j = 0
    execute_on = TIMESTEP_END
  []
  [stress_yz]
    type = ADRankTwoAux
    variable = stress_yz
    rank_two_tensor = global_stress_t_points_0 
    index_i = 1
    index_j = 2
    execute_on = TIMESTEP_END
  []
  [strain_yz]
    type = ADRankTwoAux
    variable = strain_yz
    rank_two_tensor = t_points_0_total_strain
    index_i = 1
    index_j = 2
    execute_on = TIMESTEP_END
  []
  [stress_zy]
    type = ADRankTwoAux
    variable = stress_zy
    rank_two_tensor = global_stress_t_points_0
    index_i = 2
    index_j = 1
    execute_on = TIMESTEP_END
  []
  [strain_zy]
    type = ADRankTwoAux
    variable = strain_zy
    rank_two_tensor = t_points_0_total_strain
    index_i = 2
    index_j = 1
    execute_on = TIMESTEP_END
  []
[]

[Preconditioning]
  # [./smp]
  #   type = SMP
  #   full = true
  # [../]
  [./FDP_jfnk]
    type = FDP
  [../]
[]

[Executioner]
  type = Transient
  solve_type = FD
  line_search = 'none'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  # petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount'
  # petsc_options_value = 'lu NONZERO   1e1'
  petsc_options = '-ksp_view_pmat'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-8
  dt = 1.0
  dtmin = 1.0
  end_time = 1.0
[]

[Kernels]
  [solid_disp_x]
    type = ADStressDivergenceShell2
    component = 0
    variable = disp_x
    save_in = react_disp_x
    through_thickness_order = SECOND
  []
  [solid_disp_y]
    type = ADStressDivergenceShell2
    component = 1
    variable = disp_y
    save_in = react_disp_y
    through_thickness_order = SECOND
  []
  [solid_disp_z]
    type = ADStressDivergenceShell2
    component = 2
    variable = disp_z
    save_in = react_disp_z
    through_thickness_order = SECOND
  []
  [solid_rot_x]
    type = ADStressDivergenceShell2
    component = 3
    variable = rot_x
    save_in = react_rot_x
    through_thickness_order = SECOND
    penalty = 0
  []
  [solid_rot_y]
    type = ADStressDivergenceShell2
    component = 4
    variable = rot_y
    save_in = react_rot_y
    through_thickness_order = SECOND
    penalty = 0
  []
  [solid_rot_z]
    type = ADStressDivergenceShell2
    component = 5
    variable = rot_z
    save_in = react_rot_z
    through_thickness_order = SECOND
    penalty = 0
  []
[]

[Materials]
  [elasticity_shell]
    type = ADComputeIsotropicElasticityTensorShell
    youngs_modulus = 1e6
    poissons_ratio = 0.0
    through_thickness_order = SECOND
    block = 0
  []
  [strain_shell]
    type = ADComputeIncrementalShellStrain2
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    thickness = 0.01
    through_thickness_order = SECOND
    block = 0
  []
  [stress_shell]
    type = ADComputeShellStress2
    through_thickness_order = SECOND
    block = 0
    # outputs = exodus
  []
[]

[Postprocessors]
  [xdisp_1]
    type = PointValue
    point = '1 0 0'
    variable = disp_x
  []
  [xdisp_2]
    type = PointValue
    point = '1 1 0'
    variable = disp_x
  []
  [xreact_left]
    type = NodalSum
    boundary = '3'
    variable = react_disp_x
  []
  [xreact_right]
    type = NodalSum
    boundary = '1'
    variable = react_disp_x
  []

  [stress_xx]
    type = ElementalVariableValue
    variable = 'stress_xx'
    elementid = 0
  []
  # [stress_yy]
  #   type = ElementalVariableValue
  #   variable = 'stress_yy'
  #   elementid = 0
  # []
  # [stress_xy]
  #   type = ElementalVariableValue
  #   variable = 'stress_xy'
  #   elementid = 0
  # []
  [strain_xx]
    type = ElementalVariableValue
    variable = 'strain_xx'
    elementid = 0
  []
[]

[Outputs]
  exodus = true
[]