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


int deps::require(int feature_id, bool silent /* default: false */) {  // Enable given feature
  // fails silently if requested by flag
  // silent flag is passed recursively to deps of this feature

  int res, i, j;
  bool ok;
  feature *f = features()[feature_id];
  feature_state *fs = feature_states[feature_id];

  cvm::log(description + " requiring " + f->description);

  if (!fs->available) {
    if (!silent) {
      cvm::error("Feature unavailable: \"" + f->description + "\" in object \"" + description + "\"");
    }
    return COLVARS_ERROR;
  }

  if (fs->enabled) {
    // Do not try to solve deps if already enabled
    return COLVARS_OK;
  }

  // 1) solve internal deps (self)
  for (i=0; i<f->requires_self.size(); i++) {
    res = require(f->requires_self[i], silent);
    if (res != COLVARS_OK) {
      if (!silent) {
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
      res = require(g, true);  // fail silently
      if (res == COLVARS_OK) {
        ok = true;
        break;
      }
    }
    if (!ok) {
      if (!silent) {
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
      res = children[j]->require(g, silent);
      if (res != COLVARS_OK) {
        if (!silent) {
          cvm::error("Unsatisfied dependency for \"" + f->description + "\" in object \"" + description + "\"");
        }
        return res;
      }
    }
  }

  // Actually enable feature only once everything checks out
  fs->enabled = true;
  return COLVARS_OK;
}


//     disable() {
//
//       // we need refs to parents to walk up the deps tree!
//       // or refresh
//     }
