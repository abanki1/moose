# [Mesh]
#   [mesh]
#     type = FileMeshGenerator
#     file = flatplates_xy.e
#   []
# []
[Mesh]
  [gmg]
    type = GeneratedMeshGenerator #In 2D, bottom =0, right = 1, top = 2, left = 3
    dim = 2
    nx = 1
    ny = 1
    xmin = 0.0
    xmax = 1.0
    ymin = 0.0
    ymax = 1.0
    show_info = true
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
  [stress_yz]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_xz]
    order = CONSTANT
    family = MONOMIAL
  []
  [strain_xx]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[BCs]
  [xy_fix_x]
    type = DirichletBC
    variable = disp_x
    boundary = '3' #LeftEdge
    value = 0.0
  []
  [xy_fix_y]
    type = DirichletBC
    variable = disp_y
    boundary = '0 1 2 3' #'6'#'LeftEdge'
    value = 0.0
  []
  [xy_fix_z]
    type = DirichletBC
    variable = disp_z
    boundary = '0 1 2 3' #'6' #LeftEdge
    value = 0.0
  []
  [xy_fix_rot_x]
    type = DirichletBC
    variable = rot_x
    boundary = '0 1 2 3' 
    value = 0.0
  []
  [xy_fix_rot_y]
    type = DirichletBC
    variable = rot_y
    boundary = '0 1 2 3' 
    value = 0.0
  []
  [xy_fix_rot_z]
    type = DirichletBC
    variable = rot_z
    boundary = '0 1 2 3'
    value = 0.0
  []
  # [xy_pull]
  #   type = DirichletBC
  #   variable = disp_x
  #   boundary = '1' #RightEdge
  #   value = 1.0
  # []
[]

# [DiracKernels]
#  [point1]
#    type = ConstantPointSource
#    variable = disp_x
#   #  point = '1 0 1'
#   point = '1 1 0'
#    value = 2.5 # P = 10
#  []
# []


# [AuxKernels]
#   [stress_xx]
#     type = RankTwoAux
#     rank_two_tensor = stress
#     index_i = 0
#     index_j = 0
#     variable = stress_xx
#   []
# []

[NodalKernels]
 [fx]
   type = UserForcingFunctionNodalKernel
   boundary = '1'
   function = 10
   variable = 'disp_x'
 []
[]

[Preconditioning]
  # [./smp]
  #   type = SMP
  #   full = true
  # [../]
  [FDP_jfnk]
    type = FDP
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  line_search = 'none'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  # petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount'
  # petsc_options_value = 'lu NONZERO   1e1'
  # petsc_options = '-ksp_view_pmat'
  # petsc_options = '-ksp_view_rhs'
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
  # [total_strain_xx_0]
  #   type = ADComputeIncrementalShellStrain2
  #   base_name = t_points_0
  #   rank_two_tensor = t_points_0_total_strain
  #   property_name = 'total_strain_xx_0t'
  #   displacements = 'disp_x disp_y disp_z'
  #   rotations = 'rot_x rot_y rot_z'
  #   thickness = 0.01
  #   through_thickness_order = SECOND
  #   index_i = 0
  #   index_j = 0
  #   outputs = all
  # []
  # [total_strain_xx_1]
  #   type = ADComputeIncrementalShellStrain2
  #   base_name = t_points_1
  #   rank_two_tensor = t_points_1_total_strain
  #   property_name = 'total_strain_xx_1t'
  #   displacements = 'disp_x disp_y disp_z'
  #   rotations = 'rot_x rot_y rot_z'
  #   thickness = 0.01
  #   through_thickness_order = SECOND
  #   index_i = 0
  #   index_j = 0
  #   outputs = all
  # []
#   [stress_xx_0]
#     type = ADComputeLinearElasticStress
#     rank_two_tensor = t_points_0_stress
#     property_name = 'stress_xx_0t'
#     base_name = t_points_0
#     index_i = 0
#     index_j = 0
#     outputs = all
#   []
#   [stress_xx_1]
#     type = ADComputeLinearElasticStress
#     rank_two_tensor = t_points_1_stress
#     property_name = 'stress_xx_1t'
#     base_name = t_points_1
#     index_i = 0
#     index_j = 0
#     outputs = all
#   []
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