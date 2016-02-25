#include "colvardeps.h"

/* // These are shortcuts, barely shorter than the original
void add_dep_self(int f, int dep) {
  features[f]->requires_self.push_back(dep);
}

void add_dep_alt(int f, int num, const int *deps) {
  int i;
  features[f]->requires_alt.push_back(std::vector<int>(num));
  for (i = 0; i < num; i++) {
    features[f]->requires_alt.back()[i] = deps[i];
  }
}

void add_dep_children(int) {int f, int dep) {
  features[f]->requires_children.push_back(dep);
}*/


int deps::require(int feature_id, bool dry_run /* default: false */) {  // Enable given feature
  // dry_run: fail silently, do not enable if available
  // flag is passed recursively to deps of this feature

  int res, i, j;
  bool ok;
  feature *f = features()[feature_id];
  feature_state *fs = feature_states[feature_id];

  cvm::log(description + " requiring " + f->description);

  if (fs->enabled) {
    // Do not try to solve deps if already enabled
    return COLVARS_OK;
  }

  if (!fs->available) {
    if (!dry_run) {
      cvm::error("Feature unavailable: \"" + f->description + "\" in object \"" + description + "\"");
    }
    return COLVARS_ERROR;
  }

  // 1) solve internal deps (self)
  for (i=0; i<f->requires_self.size(); i++) {
    res = require(f->requires_self[i], dry_run);
    if (res != COLVARS_OK) {
      if (!dry_run) {
          cvm::error("Unsatisfied dependency for \"" + f->description + "\" in object \"" + description + "\"");
      }
      return res;
    }
  }

  // 2) solve internal alternate deps
  for (i=0; i<f->requires_alt.size(); i++) {

    // test if one is available; if yes, enable and exit w/ success
    ok = false;
    for (j=0; j<f->requires_alt[i].size(); j++) {
      int g = f->requires_alt[i][j];
      res = require(g, true);  // see if available
      if (res == COLVARS_OK) {
        ok = true;
        if (!dry_run) require(g); // Require again, for real
        break;
      }
    }
    if (!ok) {
      if (!dry_run) {
        cvm::log("In object \"" + description + "\", no dependency satisfied among alternates:");
        for (j=0; j<f->requires_alt[i].size(); j++) {
          int g = f->requires_alt[i][j];
          cvm::log(cvm::to_str(j+1) + ". " + features()[g]->description);
        }
      }
      return COLVARS_ERROR;
    }
  }

  // 3) solve deps in children
  for (i=0; i<f->requires_children.size(); i++) {
    int g = f->requires_children[i];
//     cvm::log("children " + cvm::to_str(g));
    for (j=0; j<children.size(); j++) {
//       cvm::log("child " +  children[j]->description);
      cvm::increase_depth();
      res = children[j]->require(g, dry_run);
      cvm::decrease_depth();
      if (res != COLVARS_OK) {
        if (!dry_run) {
          cvm::error("Unsatisfied dependency for \"" + f->description + "\" in object \"" + description + "\"");
        }
        return res;
      }
    }
  }

  // Actually enable feature only once everything checks out
  if (!dry_run) fs->enabled = true;
  return COLVARS_OK;
}


//     disable() {
//
//       // we need refs to parents to walk up the deps tree!
//       // or refresh
//     }
