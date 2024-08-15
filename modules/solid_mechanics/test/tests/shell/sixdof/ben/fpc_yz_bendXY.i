[Mesh]
  [mesh]
    type = FileMeshGenerator
    file = flatplates_yz.e #5-bottom 6-left 7-top 8-right
  []
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [disp_z]
  []
  [rot_x]
  []
  [rot_y]
  []
  [rot_z]
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
[]

[BCs]
    [yz_fix_x]
        type = DirichletBC
        variable = disp_x
        boundary = '6'
        value = 0.0
    []
  [yz_fix_y]
    type = DirichletBC
    variable = disp_y
    boundary = '6'
    value = 0.0
  []
  [yz_fix_z]
    type = DirichletBC
    variable = disp_z
    boundary = '6'
    value = 0.0
  []
  [yz_fix_rot_x]
    type = DirichletBC
    variable = rot_x
    boundary = '6'
    value = 0.0
  []
  [yz_fix_rot_y]
    type = DirichletBC
    variable = rot_y
    boundary = '6'
    value = 0.0
  []
  [yz_fix_rot_z]
    type = DirichletBC
    variable = rot_z
    boundary = '6'
    value = 0.0
  []
  [yz_pull_x]
    type = DirichletBC
    variable = disp_x
    boundary = '8'
    value = -0.01
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
 [fx]
   type = UserForcingFunctionNodalKernel
   boundary = '8'
   function = -3
   variable = 'disp_x'
 []
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  line_search = 'none'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
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
    penalty = 1e-6
  []
  [solid_rot_y]
    type = ADStressDivergenceShell2
    component = 4
    variable = rot_y
    save_in = react_rot_y
    through_thickness_order = SECOND
    penalty = 1e-6
  []
  [solid_rot_z]
    type = ADStressDivergenceShell2
    component = 5
    variable = rot_z
    save_in = react_rot_z
    through_thickness_order = SECOND
    penalty = 1e-6
  []
[]

[Materials]
  [elasticity_t0]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 1e6
    poissons_ratio = 0.0
    base_name = t_points_0
  []
  [elasticity_t1]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 1e6
    poissons_ratio = 0.0
    base_name = t_points_1
  []
  [strain]
    type = ADComputeIncrementalShellStrain2
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    thickness = 0.01
    through_thickness_order = SECOND
  []
  [stress_t0]
    type = ADComputeLinearElasticStress
    base_name = t_points_0
  []
  [stress_t1]
    type = ADComputeLinearElasticStress
    base_name = t_points_1
  []
  [total_strain_yy_0]
    type = ADRankTwoCartesianComponent
    rank_two_tensor = t_points_0_total_strain
    property_name = 'total_strain_yy_0t'
    index_i = 1
    index_j = 1
    outputs = all
  []
  [total_strain_yy_1]
    type = ADRankTwoCartesianComponent
    rank_two_tensor = t_points_1_total_strain
    property_name = 'total_strain_yy_1t'
    index_i = 1
    index_j = 1
    outputs = all
  []
  [stress_yy_0]
    type = ADRankTwoCartesianComponent
    rank_two_tensor = t_points_0_stress
    property_name = 'stress_yy_0t'
    index_i = 1
    index_j = 1
    outputs = all
  []
  [stress_yy_1]
    type = ADRankTwoCartesianComponent
    rank_two_tensor = t_points_1_stress
    property_name = 'stress_yy_1t'
    index_i = 1
    index_j = 1
    outputs = all
  []
[]

[Postprocessors]
  [xdisp1]
    type = PointValue
    point = '0 1 0'
    variable = disp_x
  []
  [xdisp2]
    type = PointValue
    point = '0 1 1'
    variable = disp_x
  []
  [xreact_right]
    type = NodalSum
    boundary = 8
    variable = react_disp_x
  []
  [xreact_left]
    type = NodalSum
    boundary = 6
    variable = react_disp_x
  []
#   [z_rot_react_right]
#     type = NodalSum
#     boundary = 8
#     variable = react_rot_z
#   []
#   [z_rot_react_left]
#     type = NodalSum
#     boundary = 6
#     variable = react_rot_z
#   []
[]

[Outputs]
  exodus = true
[]
