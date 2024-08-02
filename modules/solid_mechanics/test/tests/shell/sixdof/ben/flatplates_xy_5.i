[Mesh]
  [mesh]
    type = FileMeshGenerator
    file = flatplates_xy.e
  []
  [rotate]
    type = TransformGenerator
    input = mesh
    transform = ROTATE
    vector_value = '45 0 0'
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
[]

[BCs]
  [xy_fix_x]
    type = DirichletBC
    variable = disp_x
    boundary = '6'
    value = 0.0
  []
  [xy_fix_y]
    type = DirichletBC
    variable = disp_y
    boundary = '1' #BottomLeftNode
    value = 0.0
  []
  [xy_fix_z]
    type = DirichletBC
    variable = disp_z
    boundary = '6'
    value = 0.0
  []
  [xy_fix_rot_y]
    type = DirichletBC
    variable = rot_y
    boundary = '6' #LeftEdge
    value = 0.0
  []
  [xy_pull_x]
    type = DirichletBC
    variable = disp_x
    boundary = '8' #RightEdge
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

#[NodalKernels]
#  [fx]
#    type = UserForcingFunctionNodalKernel
#    boundary = '3 4'
#    function = 1
#    variable = 'disp_x'
#  []
#[]

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
    type = ADStressDivergenceShell
    component = 0
    variable = disp_x
    save_in = react_disp_x
    through_thickness_order = SECOND
  []
  [solid_disp_y]
    type = ADStressDivergenceShell
    component = 1
    variable = disp_y
    save_in = react_disp_y
    through_thickness_order = SECOND
  []
  [solid_disp_z]
    type = ADStressDivergenceShell
    component = 2
    variable = disp_z
    save_in = react_disp_z
    through_thickness_order = SECOND
  []
  [solid_rot_x]
    type = ADStressDivergenceShell
    component = 3
    variable = rot_x
    save_in = react_rot_x
    through_thickness_order = SECOND
  []
  [solid_rot_y]
    type = ADStressDivergenceShell
    component = 4
    variable = rot_y
    save_in = react_rot_y
    through_thickness_order = SECOND
  []
[]

[Materials]
  [elasticity]
    type = ADComputeIsotropicElasticityTensorShell
    youngs_modulus = 1e6
    poissons_ratio = 0.3
    through_thickness_order = SECOND
  []
  [strain]
    type = ADComputeIncrementalShellStrain
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y'
    thickness = 0.01
    through_thickness_order = SECOND
  []
  [stress]
    type = ADComputeShellStress
    through_thickness_order = SECOND
  []
  [stress_xx_0]
    type = RankTwoCartesianComponent
    rank_two_tensor = global_stress_t_points_0
    property_name = 'stress_xx_0t'
    index_i = 0
    index_j = 0
    outputs = all
  []
  [stress_xx_1]
    type = RankTwoCartesianComponent
    rank_two_tensor = global_stress_t_points_1
    property_name = 'stress_xx_1t'
    index_i = 0
    index_j = 0
    outputs = all
  []
[]

[Postprocessors]
  [xdisp1]
    type = PointValue
    point = '1 0 0'
    variable = disp_x
  []
  [xdisp2]
    type = PointValue
    point = '1 1 0'
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
[]

[Outputs]
  exodus = true
[]
