# Test for displacement of pinched cylinder with load in Y direction
# Ref: Figure 10 and Table 6 from Dvorkin and Bathe, Eng. Comput., Vol. 1, 1984.

# Run with -pc_type svd -pc_svd_monitor if convergence issue

[Mesh]
  [./mesh]
    type = FileMeshGenerator
    file = cyl_1x1.e
  [../]
[]

[Variables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
  [../]
  [./rot_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./rot_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./rot_z]
    order = FIRST
    family = LAGRANGE
  [../]
[]
# [ICs]
#     [disp_x]
#       type = RandomIC
#       variable = disp_x
#       min = -0.01
#       max = 0.01
#     []
#     [disp_y]
#       type = RandomIC
#       variable = disp_y
#       min = -0.01
#       max = 0.01
#     []
#     [disp_z]
#       type = RandomIC
#       variable = disp_z
#       min = -0.01
#       max = 0.01
#     []
  
#     [rot_x]
#       type = RandomIC
#       variable = rot_x
#       min = -0.01
#       max = 0.01
#     []
#     [rot_y]
#       type = RandomIC
#       variable = rot_y
#       min = -0.01
#       max = 0.01
#     []
#     [rot_z]
#       type = RandomIC
#       variable = rot_z
#       min = -0.01
#       max = 0.01
#     []
#   []
[BCs]
  [./simply_support_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'CD AD'
    value = 0.0
  [../]
  [./simply_support_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'CD BC'
    value = 0.0
  [../]
  [./simply_support_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'CD AB'
    value = 0.0
  [../]
  [./simply_support_rot_x]
    type = DirichletBC
    variable = rot_x
    boundary = 'CD BC AB'
    value = 0.0
  [../]
  [./simply_support_rot_y]
    type = DirichletBC
    variable = rot_y
    boundary = 'CD AD AB'
    value = 0.0
  [../]
  [./simply_support_rot_z]
    type = DirichletBC
    variable = rot_z
    boundary = 'CD AD BC'
    # boundary = 'CD AD BC AB' #debugging attempts
    value = 0.0
  [../]
[]

[NodalKernels]
  [pinch]
    type = UserForcingFunctionNodalKernel
    boundary = 'AD' #'11'
    function = -2.5
    variable = disp_y
  []
  [constraint]
    type = PenaltyDirichletNodalKernel
    variable = rot_z
    value = 0
    penalty = 1e6
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
  petsc_options = '-ksp_view_pmat'
  nl_rel_tol = 1e-8 #was 1e-10 previously
#   nl_abs_tol = 1e-8
  dt = 1.0
  dtmin = 1.0
  end_time = 1.0
[]

[Kernels]
  [./solid_disp_x]
    type = ADStressDivergenceShell2
    block = '100'
    component = 0
    variable = disp_x
    through_thickness_order = SECOND
  [../]
  [./solid_disp_y]
    type = ADStressDivergenceShell2
    block = '100'
    component = 1
    variable = disp_y
    through_thickness_order = SECOND
  [../]
  [./solid_disp_z]
    type = ADStressDivergenceShell2
    block = '100'
    component = 2
    variable = disp_z
    through_thickness_order = SECOND
  [../]
  [./solid_rot_x]
    type = ADStressDivergenceShell2
    block = '100'
    component = 3
    variable = rot_x
    through_thickness_order = SECOND
    penalty = 0
  [../]
  [./solid_rot_y]
    type = ADStressDivergenceShell2
    block = '100'
    component = 4
    variable = rot_y
    through_thickness_order = SECOND
    penalty = 0
  [../]
  [./solid_rot_z]
    type = ADStressDivergenceShell2
    block = '100'
    component = 5
    variable = rot_z
    through_thickness_order = SECOND
    penalty = 1e6
  [../]
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
  [./strain]
    type = ADComputeIncrementalShellStrain2
    block = '100'
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
[]

[Postprocessors]
    [disp_x]
      type = PointValue
      point = '1 0 1'
      variable = disp_x
    []
    [disp_y]
      type = PointValue
      point = '0 1 1'
      variable = disp_y
    []
  []
[Postprocessors]
  [./disp_x]
    type = SideAverageValue
    boundary = 'BC'
    variable = disp_x
    
  [../]
  [./disp_y]
    type = SideAverageValue
    boundary = 'AD'
    variable = disp_y
  [../]
[]

[Outputs]
  exodus = true
[]