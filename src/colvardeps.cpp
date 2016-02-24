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


int deps::enable(int feature_id, bool silent /* default: false */) {  // Enable given feature 
  // fails silently if requested by flag
  // silent flag is passed recursively to deps of this feature
  
  int res, i, j;
  feature &f = features[feature_id];

  if (!f->available) {
    if (!silent) {
      cvm::error("Feature unavailable: " + f->description + " in object " + description);
    }
    return COLVARS_ERROR;
  }

  if (f->enabled) {
    // Do not try to solve deps if already enabled
    return COLVARS_OK;
  }

  // 1) solve internal deps (self)
  for (i=0; i<f->requires_self.size(); i++) {
    res = enable(f->requires_self[i], silent);
    if (res != COLVARS_OK) {
      return res;
    }
  }

  // 2) solve internal alternate deps
  for (i=0; i<f->requires_alt.size(); i++) {
    
    // test if one is available; if yes, enable and exit w/ success
    ok = false;
    for (j=0; j<f->requires_alt[i].size(); j++) {
      int g = f->requires_alt[i][j];
      res = enable(g, true);  // fail silently
      if (res == COLVARS_OK) {
        ok = true;
        break;
      }
    }
    if (!ok) {
      if (!silent) {
        cvm::log("No dependency satisfied among alternates:");
        for (j=0; j<f->requires_alt[i].size(); j++) {
          int g = f->requires_alt[i][j];
          cvm::log(cvm::to_str(j) + ". " + features[g]->description);
        }
      }
      return COLVARS_ERROR;
    }
  }

  // 3) solve deps in children
  for (i=0; i<f->requires_children.size(); i++) {
    int g = f->requires_alt[i];
    for (j=0; j<children.size(); j++) {
      res = children[j].enable(g, silent);
      if (res != COLVARS_OK) {
        if (!silent) {
          cvm::error("Unsatisfied dependency for " + f->description + " in object " + description);
        }
        return res;
      }
    }
  }

  // Actually enable feature only once everything checks out
  f->enabled = true;
  return COLVARS_OK;
}


//     disable() {
// 
//       // we need refs to parents to walk up the deps tree!
//       // or refresh
//     }
