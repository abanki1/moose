//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PNSFVSolidHeatTransferPhysics.h"
#include "WCNSFVCoupledAdvectionPhysicsHelper.h"
#include "WCNSFVFlowPhysics.h"
#include "NSFVAction.h"

registerPhysicsBaseTasks("NavierStokesApp", PNSFVSolidHeatTransferPhysics);
registerMooseAction("NavierStokesApp", PNSFVSolidHeatTransferPhysics, "add_variable");
registerMooseAction("NavierStokesApp", PNSFVSolidHeatTransferPhysics, "add_ic");
registerMooseAction("NavierStokesApp", PNSFVSolidHeatTransferPhysics, "add_fv_kernel");
registerMooseAction("NavierStokesApp", PNSFVSolidHeatTransferPhysics, "add_fv_bc");
registerMooseAction("NavierStokesApp", PNSFVSolidHeatTransferPhysics, "add_material");
registerMooseAction("NavierStokesApp", PNSFVSolidHeatTransferPhysics, "add_preconditioning");

InputParameters
PNSFVSolidHeatTransferPhysics::validParams()
{
  InputParameters params = HeatConductionFV::validParams();
  params.addClassDescription("Define the Navier Stokes porous media solid energy equation");

  // Swap out some parameters, base class is not specific to porous media
  // Variables
  params.renameParam("temperature_name",
                     "solid_temperature_variable",
                     "Name of the solid phase temperature variable");
  params.set<VariableName>("solid_temperature_variable") = NS::T_solid;
  params.addParam<NonlinearVariableName>(
      "fluid_temperature_variable", NS::T_fluid, "Name of the fluid temperature variable");
  MooseEnum face_interpol_types("average skewness-corrected", "average");
  params.addParam<MooseEnum>(
      "solid_temperature_face_interpolation",
      face_interpol_types,
      "The numerical scheme to interpolate the temperature/energy to the "
      "face for conduction (separate from the advected quantity interpolation).");
  params.addParam<bool>(
      "solid_temperature_two_term_bc_expansion",
      true,
      "If a two-term Taylor expansion is needed for the determination of the boundary values"
      "of the temperature/energy.");

  // Porous media parameters
  // TODO: ensure consistency with fluid energy physics
  params.transferParam<MooseFunctorName>(NSFVAction::validParams(), "porosity");

  // Material properties
  params.suppressParameter<MaterialPropertyName>("specific_heat");
  params.addParam<MooseFunctorName>("cp_solid", NS::cp + "_solid", "Specific heat functor");
  params.suppressParameter<MaterialPropertyName>("density");
  params.addParam<MooseFunctorName>("rho_solid", NS::density + "_solid", "Density functor");
  params.addParam<std::vector<std::vector<SubdomainName>>>(
      "thermal_conductivity_blocks", "Blocks which each thermal conductivity is defined");
  params.suppressParameter<MooseFunctorName>("thermal_conductivity_functor");
  params.addRequiredParam<std::vector<MooseFunctorName>>(
      "thermal_conductivity_solid",
      "Thermal conductivity, which may have different names depending on the subdomain");

  // Ambient convection with the liquid phase parameters
  params.addParam<std::vector<std::vector<SubdomainName>>>(
      "ambient_convection_blocks", {}, "The blocks where the ambient convection is present.");
  params.addParam<std::vector<MooseFunctorName>>(
      "ambient_convection_alpha",
      {},
      "The heat exchange coefficients for each block in 'ambient_convection_blocks'.");
  params.addParam<std::vector<MooseFunctorName>>(
      "ambient_convection_temperature",
      {NS::T_fluid},
      "The fluid temperature for each block in 'ambient_convection_blocks'.");

  // Heat source in solid porous medium parameters
  params.addParam<std::vector<SubdomainName>>("external_heat_source_blocks",
                                              std::vector<SubdomainName>(),
                                              "The blocks where the heat source is present.");
  params.addParam<MooseFunctorName>(
      "external_heat_source",
      "The name of a functor which contains the external heat source for the energy equation.");
  params.addParam<Real>(
      "external_heat_source_coeff", 1.0, "Multiplier for the coupled heat source term.");
  params.addParam<bool>("use_external_enthalpy_material",
                        false,
                        "To indicate if the enthalpy material is set up outside of the action.");

  // Numerical scheme
  params.addParam<unsigned short>(
      "ghost_layers", 2, "Number of layers of elements to ghost near process domain boundaries");
  // Preconditioning has not been derived for NSFV + porous heat transfer at this point
  MooseEnum pc_options("default none", "none");
  params.set<MooseEnum>("preconditioning") = pc_options;
  params.suppressParameter<MooseEnum>("preconditioning");

  // Parameter groups
  params.addParamNamesToGroup("rho_solid cp_solid thermal_conductivity_solid "
                              "thermal_conductivity_blocks use_external_enthalpy_material",
                              "Material properties");
  params.addParamNamesToGroup("ambient_convection_alpha ambient_convection_blocks "
                              "ambient_convection_temperature",
                              "Ambient convection");
  params.addParamNamesToGroup(
      "external_heat_source_blocks external_heat_source external_heat_source_coeff",
      "Solid porous medium heat source");
  params.addParamNamesToGroup(
      "solid_temperature_face_interpolation solid_temperature_two_term_bc_expansion",
      "Numerical scheme");
  params.addParamNamesToGroup("ghost_layers", "Advanced");

  return params;
}

PNSFVSolidHeatTransferPhysics::PNSFVSolidHeatTransferPhysics(const InputParameters & parameters)
  : HeatConductionFV(parameters),
    _solid_temperature_name(getParam<VariableName>("solid_temperature_variable")),
    _fluid_temperature_name(getParam<NonlinearVariableName>("fluid_temperature_variable")),
    _porosity_name(getParam<MooseFunctorName>(NS::porosity)),
    _density_name(getParam<MooseFunctorName>("rho_solid")),
    _specific_heat_name(getParam<MooseFunctorName>("cp_solid")),
    _thermal_conductivity_blocks(
        parameters.isParamValid("thermal_conductivity_blocks")
            ? getParam<std::vector<std::vector<SubdomainName>>>("thermal_conductivity_blocks")
            : std::vector<std::vector<SubdomainName>>()),
    _thermal_conductivity_name(
        getParam<std::vector<MooseFunctorName>>("thermal_conductivity_solid")),
    _ambient_convection_blocks(
        getParam<std::vector<std::vector<SubdomainName>>>("ambient_convection_blocks")),
    _ambient_convection_alpha(getParam<std::vector<MooseFunctorName>>("ambient_convection_alpha")),
    _ambient_temperature(getParam<std::vector<MooseFunctorName>>("ambient_convection_temperature"))
{
  saveNonlinearVariableName(_solid_temperature_name);

  // Parameter checks
  if (getParam<std::vector<MooseFunctorName>>("ambient_convection_temperature").size() != 1)
    checkVectorParamsSameLengthIfSet<MooseFunctorName, MooseFunctorName>(
        "ambient_convection_alpha", "ambient_convection_temperature");
  checkSecondParamSetOnlyIfFirstOneSet("external_heat_source", "external_heat_source_coeff");
}

void
PNSFVSolidHeatTransferPhysics::addNonlinearVariables()
{
  // Dont add if the user already defined the variable
  if (nonlinearVariableExists(_solid_temperature_name,
                              /*error_if_aux=*/true))
    checkBlockRestrictionIdentical(_solid_temperature_name,
                                   getProblem().getVariable(0, _solid_temperature_name).blocks());
  else
  {
    auto params = getFactory().getValidParams("INSFVEnergyVariable");
    assignBlocks(params, _blocks);
    params.set<std::vector<Real>>("scaling") = {getParam<Real>("temperature_scaling")};
    params.set<MooseEnum>("face_interp_method") =
        getParam<MooseEnum>("solid_temperature_face_interpolation");
    params.set<bool>("two_term_boundary_expansion") =
        getParam<bool>("solid_temperature_two_term_bc_expansion");
    getProblem().addVariable("INSFVEnergyVariable", _solid_temperature_name, params);
  }
}

void
PNSFVSolidHeatTransferPhysics::addFVKernels()
{
  if (isTransient())
    addPINSSolidEnergyTimeKernels();

  addPINSSolidEnergyHeatConductionKernels();
  if (getParam<std::vector<MooseFunctorName>>("ambient_convection_alpha").size())
    addPINSSolidEnergyAmbientConvection();
  if (isParamValid("external_heat_source"))
    addPINSSolidEnergyExternalHeatSource();
}

void
PNSFVSolidHeatTransferPhysics::addPINSSolidEnergyTimeKernels()
{
  const auto kernel_type = "PINSFVEnergyTimeDerivative";
  const auto kernel_name = prefix() + "pins_solid_energy_time";

  InputParameters params = getFactory().getValidParams(kernel_type);
  assignBlocks(params, _blocks);
  params.set<NonlinearVariableName>("variable") = _solid_temperature_name;
  params.set<MooseFunctorName>(NS::density) = _density_name;

  // The '_solid' suffix has been declared when creating the INSFVEnthalpyMaterial
  // only for thermal functor material properties
  // Using this derivative we can model non-constant specific heat
  if (getProblem().hasFunctor(NS::time_deriv(NS::specific_enthalpy) + "_solid",
                              /*thread_id=*/0))
    params.set<MooseFunctorName>(NS::time_deriv(NS::specific_enthalpy)) =
        NS::time_deriv(NS::specific_enthalpy) + "_solid";
  else
    params.set<MooseFunctorName>(NS::cp) = _specific_heat_name;

  params.set<MooseFunctorName>(NS::porosity) = _porosity_name;
  // If modeling a variable density
  if (getProblem().hasFunctor(NS::time_deriv(_density_name),
                              /*thread_id=*/0))
  {
    params.set<MooseFunctorName>(NS::time_deriv(NS::density)) = NS::time_deriv(_density_name);
    params.set<MooseFunctorName>(NS::specific_enthalpy) = NS::specific_enthalpy + "_solid";
  }
  params.set<bool>("is_solid") = true;

  getProblem().addFVKernel(kernel_type, kernel_name, params);
}

void
PNSFVSolidHeatTransferPhysics::addPINSSolidEnergyHeatConductionKernels()
{
  const auto vector_conductivity = processThermalConductivity();

  const auto kernel_type =
      vector_conductivity ? "PINSFVEnergyAnisotropicDiffusion" : "PINSFVEnergyDiffusion";

  InputParameters params = getFactory().getValidParams(kernel_type);
  params.set<NonlinearVariableName>("variable") = _solid_temperature_name;
  params.set<MooseFunctorName>(NS::porosity) = _porosity_name;

  // Set block restrictions
  const bool combined = _thermal_conductivity_blocks.size() > 1;
  std::vector<SubdomainName> thermal_conductivity_blocks;
  for (const auto & block_group : _thermal_conductivity_blocks)
    thermal_conductivity_blocks.insert(thermal_conductivity_blocks.end(),
                                       std::make_move_iterator(block_group.begin()),
                                       std::make_move_iterator(block_group.end()));
  const auto block_names =
      _thermal_conductivity_blocks.size() ? thermal_conductivity_blocks : _blocks;
  assignBlocks(params, block_names);

  // Set thermal conductivity
  const auto conductivity_name = vector_conductivity ? NS::kappa : NS::k;
  if (combined)
    params.set<MooseFunctorName>(conductivity_name) = prefix() + "combined_thermal_conductivity";
  else
    params.set<MooseFunctorName>(conductivity_name) = _thermal_conductivity_name[0];

  getProblem().addFVKernel(kernel_type, prefix() + "pins_energy_diffusion", params);
}

void
PNSFVSolidHeatTransferPhysics::addPINSSolidEnergyAmbientConvection()
{
  const auto num_convection_blocks = _ambient_convection_blocks.size();
  const auto num_used_blocks = num_convection_blocks ? num_convection_blocks : 1;

  const auto kernel_type = "PINSFVEnergyAmbientConvection";
  InputParameters params = getFactory().getValidParams(kernel_type);
  params.set<NonlinearVariableName>("variable") = _solid_temperature_name;
  params.set<MooseFunctorName>(NS::T_solid) = _solid_temperature_name;
  params.set<bool>("is_solid") = true;

  for (const auto block_i : make_range(num_used_blocks))
  {
    std::string block_name = "";
    if (num_convection_blocks)
    {
      params.set<std::vector<SubdomainName>>("block") = _ambient_convection_blocks[block_i];
      block_name = Moose::stringify(_ambient_convection_blocks[block_i]);
    }
    else
    {
      assignBlocks(params, _blocks);
      block_name = std::to_string(block_i);
    }

    params.set<MooseFunctorName>("h_solid_fluid") = _ambient_convection_alpha[block_i];
    if (_ambient_temperature.size() > 1)
      params.set<MooseFunctorName>(NS::T_fluid) = _ambient_temperature[block_i];
    else
      params.set<MooseFunctorName>(NS::T_fluid) = _ambient_temperature[0];

    getProblem().addFVKernel(kernel_type, prefix() + "ambient_convection_" + block_name, params);
  }
}

void
PNSFVSolidHeatTransferPhysics::addPINSSolidEnergyExternalHeatSource()
{
  const std::string kernel_type = "FVCoupledForce";
  InputParameters params = getFactory().getValidParams(kernel_type);
  params.set<NonlinearVariableName>("variable") = _solid_temperature_name;
  const auto & source_blocks = getParam<std::vector<SubdomainName>>("external_heat_source_blocks");
  if (source_blocks.size())
    assignBlocks(params, source_blocks);
  else
    assignBlocks(params, _blocks);
  params.set<MooseFunctorName>("v") = getParam<MooseFunctorName>("external_heat_source");
  params.set<Real>("coef") = getParam<Real>("external_heat_source_coeff");

  getProblem().addFVKernel(kernel_type, prefix() + "external_heat_source", params);
}

bool
PNSFVSolidHeatTransferPhysics::processThermalConductivity()
{
  if (isParamValid("thermal_conductivity_blocks"))
    checkBlockwiseConsistency<MooseFunctorName>("thermal_conductivity_blocks",
                                                {"thermal_conductivity_solid"});
  bool have_scalar = false;
  bool have_vector = false;

  for (const auto i : index_range(_thermal_conductivity_name))
  {
    // First, check if the name is just a number (only in case of isotropic conduction)
    if (MooseUtils::parsesToReal(_thermal_conductivity_name[i]))
      have_scalar = true;
    // Now we determine what kind of functor we are dealing with
    else
    {
      if (getProblem().hasFunctorWithType<ADReal>(_thermal_conductivity_name[i],
                                                  /*thread_id=*/0))
        have_scalar = true;
      else
      {
        if (getProblem().hasFunctorWithType<ADRealVectorValue>(_thermal_conductivity_name[i],
                                                               /*thread_id=*/0))
          have_vector = true;
        else
          paramError("thermal_conductivity_solid",
                     "We only allow functor of type (AD)Real or (AD)RealVectorValue for thermal "
                     "conductivity! Functor '" +
                         _thermal_conductivity_name[i] + "' is not of the requested type.");
      }
    }
  }

  if (have_vector == have_scalar)
    paramError("thermal_conductivity_solid",
               "The entries on thermal conductivity shall either be scalars or vectors, mixing "
               "them is not supported!");
  return have_vector;
}

void
PNSFVSolidHeatTransferPhysics::addMaterials()
{
  if (!getParam<bool>("use_external_enthalpy_material"))
  {
    InputParameters params = getFactory().getValidParams("INSFVEnthalpyFunctorMaterial");
    assignBlocks(params, _blocks);

    params.set<MooseFunctorName>(NS::density) = _density_name;
    params.set<MooseFunctorName>(NS::cp) = _specific_heat_name;
    params.set<MooseFunctorName>("temperature") = _solid_temperature_name;
    params.set<MaterialPropertyName>("declare_suffix") = "solid";

    getProblem().addMaterial(
        "INSFVEnthalpyFunctorMaterial", prefix() + "ins_enthalpy_material", params);
  }

  // Combine the functors (combining scalars and vectors is not currently supported)
  if (_thermal_conductivity_name.size() > 1)
  {
    const auto vector_conductivity = processThermalConductivity();
    const auto combiner_functor = vector_conductivity ? "PiecewiseByBlockVectorFunctorMaterial"
                                                      : "PiecewiseByBlockFunctorMaterial";
    InputParameters params = getFactory().getValidParams(combiner_functor);
    params.set<MooseFunctorName>("prop_name") = prefix() + "combined_thermal_conductivity";
    std::vector<SubdomainName> blocks_list;
    std::map<std::string, std::string> blocks_to_functors;
    for (const auto i : index_range(_thermal_conductivity_name))
    {
      for (const auto & block : _thermal_conductivity_blocks[i])
      {
        blocks_list.push_back(block);
        blocks_to_functors.insert(
            std::pair<std::string, std::string>(block, _thermal_conductivity_name[i]));
      }
    }
    params.set<std::vector<SubdomainName>>("block") = blocks_list;
    params.set<std::map<std::string, std::string>>("subdomain_to_prop_value") = blocks_to_functors;
    getProblem().addMaterial(combiner_functor, prefix() + "thermal_conductivity_combiner", params);
  }
}

InputParameters
PNSFVSolidHeatTransferPhysics::getAdditionalRMParams() const
{
  unsigned short necessary_layers = getParam<unsigned short>("ghost_layers");
  if (getParam<MooseEnum>("solid_temperature_face_interpolation") == "skewness-corrected")
    necessary_layers = std::max(necessary_layers, (unsigned short)3);

  // Just an object that has a ghost_layers parameter and performs geometric, algebraic, and
  // coupling ghosting
  const std::string kernel_type = "INSFVMixingLengthReynoldsStress";
  InputParameters params = getFactory().getValidParams(kernel_type);
  params.template set<unsigned short>("ghost_layers") = necessary_layers;

  return params;
}
