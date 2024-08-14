# Test for displacement of pinched cylinder
# Ref: Figure 10 and Table 6 from Dvorkin and Bathe, Eng. Comput., Vol. 1, 1984.

# A cylinder of radius 1 m and length 2 m (along Z axis) with clamped ends
# (at z = 0 and 2 m) is pinched at mid-length by placing point loads of 10 N
# at (1, 0, 1) and (-1, 0, 1). Due to the symmetry of the problem, only 1/8th
# of the cylinder needs to be modeled.

# The normalized series solution for the displacement at the loading point is
# w = Wc E t / P = 164.24; where Wc is the displacement in m, E is the Young's
# modulus, t is the thickness and P is the point load.
# w = 164.24 * 1e6 * 0.01 / 2.5 =

# For this problem, E = 1e6 Pa, L = 2 m, R = 1 m, t = 0.01 m, P = 10 N and
# Poisson's ratio = 0.3. FEM results from different mesh discretizations are
# presented below. Only the 10x10 mesh is included as a test.

# Mesh of 1/8 cylinder |  FEM/analytical (Moose) | FEM/analytical (Dvorkin)
#                      |ratio of normalized disp.| ratio of normalized disp.
#----------------------|-------------------------|-------------------------
#     10 x 10          |          0.806          |        0.83
#     20 x 20          |          1.06           |        0.96
#     40 x 40          |          0.95           |         -
#     80 x 160         |          0.96           |         -

# The results from FEM analysis matches well with the series solution and with
# the solution presented by Dvorkin and Bathe (1984).

# Run with -pc_type svd -pc_svd_monitor if convergence issue

# [GlobalParams]
#   use_displaced_mesh = true
# []

[Mesh]
  [mesh]
    type = FileMeshGenerator
    file = cyl_1x2.e
  []
  [all_nodes]
    type = BoundingBoxNodeSetGenerator
    input = mesh
    top_right = '1e6 1e6 1e6'
    bottom_left = '-1e6 -1e6 -1e6'
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

[BCs]
  [simply_support_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'CD AD'
    # boundary ='CD'
    value = 0.0
  []
  [simply_support_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'CD BC'
    # boundary ='CD'
    value = 0.0
  []
  [simply_support_z]
    type = DirichletBC
    variable = disp_z
    boundary ='CD AB'
    value = 0.0
  []
  [simply_support_rot_x]
    type = DirichletBC
    variable = rot_x
    # boundary = 'CD AB BC'
    boundary ='CD BC'
    # boundary = all_nodes
    value = 0.0
  []
  [simply_support_rot_y]
    type = DirichletBC
    variable = rot_y
    # boundary = 'CD AB AD'
    # boundary ='CD AD'
    boundary = all_nodes
    value = 0.0
  []
  [simply_support_rot_z]
    type = DirichletBC
    variable = rot_z
    # boundary = 'CD AD BC'
    # boundary ='CD'
    boundary = all_nodes
    value = 0.0
  []
[]

[NodalKernels]
  [pinch]
    type = UserForcingFunctionNodalKernel
    boundary = 'BC' #'10'
    function = -2.5
    variable = disp_x
  []
  # [./constraint_x]
  #   type = PenaltyDirichletNodalKernel
  #   variable = rot_x
  #   value = 0
  #   penalty = 1e12
  # []
  # [./constraint_y]
  #   type = PenaltyDirichletNodalKernel
  #   variable = rot_y
  #   value = 0
  #   penalty = 1e12
  # []
  [./constraint_z]
    type = PenaltyDirichletNodalKernel
    variable = rot_z
    value = 0
    penalty = 1e12
  []
[]

[Preconditioning]
#   [./smp]
#     type = SMP
#     full = true
#   [../]
  [FDP_jfnk]
    type = FDP
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  # line_search = 'none'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  petsc_options = '-ksp_view_pmat'
  nl_rel_tol = 1e-8
#   nl_abs_tol = 1e-8
  dt = 1.0
  dtmin = 1.0
  end_time = 1.0
[]

[Kernels]
  [solid_disp_x]
    type = ADStressDivergenceShell2
    block = '100'
    component = 0
    variable = disp_x
    through_thickness_order = SECOND
  []
  [solid_disp_y]
    type = ADStressDivergenceShell2
    block = '100'
    component = 1
    variable = disp_y
    through_thickness_order = SECOND
  []
  [solid_disp_z]
    type = ADStressDivergenceShell2
    block = '100'
    component = 2
    variable = disp_z
    through_thickness_order = SECOND
  []
  [solid_rot_x]
    type = ADStressDivergenceShell2
    block = '100'
    component = 3
    variable = rot_x
    through_thickness_order = SECOND
    penalty = 1e6
  []
  [solid_rot_y]
    type = ADStressDivergenceShell2
    block = '100'
    component = 4
    variable = rot_y
    through_thickness_order = SECOND
    penalty = 1e6
  []
  [solid_rot_z]
    type = ADStressDivergenceShell2
    block = '100'
    component = 5
    variable = rot_z
    through_thickness_order = SECOND
    penalty = 1e6
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

[Outputs]
  exodus = true
[]