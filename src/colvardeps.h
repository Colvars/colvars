#include "colvarmodule.h"

#ifndef COLVARDEPS_H
#define COLVARDEPS_H

/// Parent class for a member object of a bias, cv or cvc etc. containing dependencies
/// (features) and handling dependency resolution


// OPTION A:
// Simpler solution: update the "features" member object in the {colvarbias, colvar, cvc, atomgroup} constructor/initializer
// No need for derived classes, code goes into colvar etc. classes
// features set in cvc constructor, then refined in colvar_angle, etc.


// OPTION B:
// To be implemented as colvar_deps, cvc_deps, atomgroup_deps ?
// deps ->    colvar_deps -> colvar
//   \ \->       cvc_deps -> cvc
//    \--> atomgroup_deps -> atomgroup



// Some features like colvar::f_linear have no dependencies, require() doesn't enable anything but fails if unavailable
// Policy: those features should be either enabled or unavailable at all times

// It seems important to have available default to false (for safety) and enabled to false (for efficiency)

// Features need to be in a common namespace so that parents can list the features of children
// derived class features can be added beyond the range of base class features (eg. 1site_force), to be used for self deps

class deps {
public:

  deps() {}
  ~deps() {
    // Do not delete features if it's static
//     for (i=0; i<features.size(); i++) {
//       if (features[i] != NULL) delete features[i];
//     }
  }

  // Subclasses should initialize the following members:

  std::string description; // reference to object name (cv, cvc etc.)

  /// This contains the current state of each feature for each object
  struct feature_state {
    /// Available means supported, subject to dependencies as listed
    /// unavailable until declared otherwise
    bool available = false;
    /// Currently enabled - this flag is subject to change dynamically
    /// TODO consider implications for dependency solving: anyone who disables
    /// it should trigger a refresh of parent objects
    bool enabled = false;     // see if this should be private depending on implementation
    // bool enabledOnce; // this should trigger an update when object is evaluated
  };

  /// List of the state of all features
  std::vector<feature_state *> feature_states;

  /// Describes a feature and its dependecies
  /// used in a static array within each subclass
  class feature {

  public:
    feature() {}
    ~feature() {}

    std::string description; // Set by derived object initializer

    // features that this feature requires in the same object
    // NOTE: we have no safety mechanism against circular dependencies, however, they would have to be internal to an object (ie. requires_self or requires_alt)
    std::vector<int> requires_self;


    // sets of features that are required in an alternate way
    // when parent feature is enabled, if none are enabled, the first one listed that is available will be enabled
    std::vector<std::vector<int>> requires_alt;

    // features that this feature requires in children
    std::vector<int> requires_children;
  };

  // Accessor to list of all features, static in each derived class
  // Can be iterated over by going through subclass-specific enums
  virtual std::vector<feature *>&features() = 0;

  // pointers to objects this object depends on
  // list should be maintained by any code that modifies the object
  // this could be secured by making lists of colvars / cvcs / atom groups private and modified through accessor functions
  std::vector<deps *> children;


  // std::vector<deps *> parents; // Needed to trigger a refresh if capabilities of this object change

  // End of members to be initialized by subclasses


  int require(int f, bool dry_run = false);  // enable a feature and recursively solve its dependencies
  // dry_run is set to true to recursively test if a feature is available, without enabling it
//     int disable(int f);

  // At this point it is unclear which of the following refresh mechanisms will be useful
  // int resolve() {} // Propagate "enables" downwards to children then trigger resolve in children
  // int refresh() {} // Trigger refresh in children, then calculate whether each feature is available; sets available flag


  /* Mockup code
  features[f_Jacobian]->available
  features[f_Jacobian]->enable();

  if (!children[i]->features[f_Jacobian]->enable()) {
    err_failed_dependency(description, features[f_Jacobian]->description, children[i]->features[f_Jacobian]->description);
    // Ideally, give parent object name, parent feature descr, and child feature that failed
  }
  */

  // Those are also documented by their description std::string
  // see initializers in derived classes
  enum features_colvar {
    /// \brief Calculate colvar value
    f_cv_value,
    /// \brief is scalar
    f_cv_scalar,
    f_cv_linear,
    f_cv_homogeneous,
    /// \brief Gradients are calculated and temporarily stored, so
    /// that external forces can be applied
    f_cv_gradient,
    /// \brief Collect atomic gradient data from all cvcs into vector
    /// atomic_gradient
    f_cv_collect_gradient,
    /// \brief Calculate the velocity with finite differences
    f_cv_fdiff_velocity,
    /// \brief The system force is calculated, projecting the atomic
    /// forces on the inverse gradient
    f_cv_system_force,
    /// \brief Calculate system force from atomic forces
    f_cv_system_force_calc,
    /// \brief The variable has a harmonic restraint around a moving
    /// center with fictitious mass; bias forces will be applied to
    /// the center
    f_cv_extended_Lagrangian,
    /// \brief The extended system coordinate undergoes Langevin dynamics
    f_cv_Langevin,
    /// \brief Output the potential and kinetic energies
    /// (for extended Lagrangian colvars only)
    f_cv_output_energy,
    /// \brief Estimate Jacobian derivative
    f_cv_Jacobian,
    /// \brief Compute Jacobian derivative explicitly
    f_cv_Jacobian_calc,
    /// \brief Report the Jacobian force as part of the system force
    /// if disabled, apply a correction internally to cancel it
    f_cv_report_Jacobian,
    /// \brief Output the value to the trajectory file (on by default)
    f_cv_output_value,
    /// \brief Output the velocity to the trajectory file
    f_cv_output_velocity,
    /// \brief Output the applied force to the trajectory file
    f_cv_output_applied_force,
    /// \brief Output the system force to the trajectory file
    f_cv_output_system_force,
    /// \brief A lower boundary is defined
    f_cv_lower_boundary,
    /// \brief An upper boundary is defined
    f_cv_upper_boundary,
    /// \brief Provide a discretization of the values of the colvar to
    /// be used by the biases or in analysis (needs lower and upper
    /// boundary)
    f_cv_grid,
    /// \brief Apply a restraining potential (|x-xb|^2) to the colvar
    /// when it goes below the lower wall
    f_cv_lower_wall,
    /// \brief Apply a restraining potential (|x-xb|^2) to the colvar
    /// when it goes above the upper wall
    f_cv_upper_wall,
    /// \brief Compute running average
    f_cv_runave,
    /// \brief Compute time correlation function
    f_cv_corrfunc,
    /// \brief Value and gradient computed by user script
    f_cv_scripted,
    /// \brief Number of colvar features
    f_cv_ntot
  };



    // Shorthand macros for describing dependencies
#define f_description(f, d) features()[f]->description = d;
#define f_req_self(f, g) features()[f]->requires_self.push_back(g);
#define f_req_children(f, g) features()[f]->requires_children.push_back(g);
#define f_req_alt2(f, g, h) features()[f]->requires_alt.push_back(std::vector<int>(2));\
  features()[f]->requires_alt.back()[0] = g;                                           \
  features()[f]->requires_alt.back()[1] = h;

  inline void init_cv_requires() {
    for (int i = 0; i < f_cv_ntot; i++) {
      features().push_back(new feature);
    }

    f_description(f_cv_value, "value")
    f_req_children(f_cv_value, f_cvc_value)

    f_description(f_cv_gradient, "gradient")
    f_req_self(f_cv_gradient, f_cv_value)
    f_req_children(f_cv_gradient, f_cvc_gradient)

    f_description(f_cv_collect_gradient, "collect gradient")
    f_req_self(f_cv_collect_gradient, f_cv_gradient)

    // System force: either trivial (spring force) through extended Lagrangian, or calculated explicitly
    f_description(f_cv_system_force, "system force")
    f_req_alt2(f_cv_system_force, f_cv_extended_Lagrangian, f_cv_system_force_calc)

    // Deps for explicit system force calculation
    f_description(f_cv_system_force_calc, "system force calculation")
    f_req_self(f_cv_system_force_calc, f_cv_scalar)
    f_req_self(f_cv_system_force_calc, f_cv_linear)
    f_req_children(f_cv_system_force_calc, f_cvc_system_force)

    // Jacobian force: either trivial (zero) through extended Lagrangian, or calculated explicitly
    f_description(f_cv_Jacobian, "Jacobian derivative")
    f_req_alt2(f_cv_Jacobian, f_cv_extended_Lagrangian, f_cv_Jacobian_calc)

    // Deps for explicit Jacobian calculation
    f_description(f_cv_Jacobian_calc, "Jacobian derivative calculation")
    f_req_self(f_cv_Jacobian_calc, f_cv_scalar)
    f_req_self(f_cv_Jacobian_calc, f_cv_linear)
    f_req_children(f_cv_Jacobian_calc, f_cvc_Jacobian)

    f_description(f_cv_report_Jacobian, "report Jacobian force")

    f_description(f_cv_output_value, "output value")
    f_req_self(f_cv_output_value, f_cv_value)

    f_description(f_cv_extended_Lagrangian, "extended Lagrangian")
    f_description(f_cv_Langevin, "Langevin dynamics")
    f_description(f_cv_linear, "linear")
    f_description(f_cv_scalar, "scalar")
  }

  enum features_cvc {
    f_cvc_value,
    f_cvc_scalar,
    f_cvc_gradient,
    f_cvc_system_force,
    f_cvc_inv_gradient,
    f_cvc_Jacobian,
    f_cvc_ntot
  };

  inline void init_cvc_requires() {
    for (int i = 0; i < deps::f_cvc_ntot; i++) {
      features().push_back(new feature);
    }

    f_description(f_cvc_value, "value")
    f_req_children(f_cvc_value, f_ag_coordinates)

    f_description(f_cvc_scalar, "scalar")

    f_description(f_cvc_gradient, "gradient")
    f_req_self(f_cvc_gradient, f_cvc_value)

    f_description(f_cvc_system_force, "system force")
    f_req_self(f_cvc_system_force, f_cvc_inv_gradient)

    f_description(f_cvc_inv_gradient, "inverse gradient")
    f_req_self(f_cvc_inv_gradient, f_cvc_gradient)

    f_description(f_cvc_Jacobian, "Jacobian")
    f_req_self(f_cvc_Jacobian, f_cvc_inv_gradient)
  }

  enum features_atomgroup {
    f_ag_coordinates,
    f_ag_fit,
    f_ag_fit_gradient_group,// TODO check that these are sometimes needed separately
                            // maybe for minimum RMSD?
    f_ag_fit_gradient_ref,
    f_ag_atom_forces,
    f_ag_ntot
  };

  inline void init_ag_requires() {
    for (int i = 0; i < f_ag_ntot; i++) {
      features().push_back(new feature);
    }

    f_description(f_ag_coordinates, "coordinates")
    f_description(f_ag_fit, "optimal fit")
    f_description(f_ag_fit_gradient_group, "fit gradient for main group")
    f_description(f_ag_fit_gradient_ref, "fit gradient for reference group")
    f_description(f_ag_atom_forces, "atomic forces")
  }


  /// \brief print all enabled features and those of children, for debugging
  void print_state() {
    cvm::log("Enabled features of object " + description);
    for (int i = 0; i<feature_states.size(); i++) {
      if (feature_states[i]->enabled)
        cvm::log(cvm::to_str(i) + " " + features()[i]->description);
    }
    for (int i=0; i<children.size(); i++) {
      cvm::log("Child " + cvm::to_str(i+1));
      cvm::increase_depth();
      children[i]->print_state();
      cvm::decrease_depth();
    }
  }
};

#endif


