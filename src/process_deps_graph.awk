#!/usr/bin/awk -f

# This script parses a Colvars C++ source files containing calls to colvardeps::init_feature()
# and print a dependency graph in dot (graphviz) format

/^ *init_feature\(/ {
  # Create a node (text box) for each feature

  if (!match($0, /init_feature\(f_([^,]*), "/, res)) exit(1)
  id = res[1]
  if (!match($0, /init_feature\(f_([^_]*)/, res)) exit(1)
  obj = res[1]
  if (!match($0, /, "([^,]*)"/, res)) exit(1)
  label = res[1]
  if (!match($0, /f_type_([^,]*)\);/, res)) 
    type = "static"
  else
    type = res[1]
    
  color["cvb"] = "red"
  color["cv"] = "black"
  color["cvc"] = "blue"
  color["ag"] = "darkgreen"

  shape["dynamic"] = "ellipse"
  shape["static"] = "box"
  shape["user"] = "box"

  style["dynamic"] = "dashed, bold"
  style["static"] = "solid"
  style["user"] = "rounded,bold"

  depstyle["self"] = ""
  depstyle["exclude"] = "[dir = \"both\" arrowhead = \"tee\" arrowtail = \"tee\" style = dotted]"
  depstyle["children"] = ""
  depstyle["alt"] = "[style = dashed]"

  print "  " id " [label = \"" label "\" fontcolor = \"" color[obj] "\" color = \"" \
    color[obj] "\" shape = \"" shape[type] "\" style = \"" style[type] "\"];"
}

/^ *require_feature_/ {
  # Create an edge for each dependency

  if (!match($0, /require_feature_([^\(]*)/, res)) exit(1)
  deptype = res[1]

  # extract feature names
  sub(/\);.*/, "")
  n = split(gensub(/[^a-z_]*f_([a-z_]*)[\);]*/, " \\1", "g"), features)

  for (i=2; i<=n; i++) {
    color[i] = color[gensub(/_.*/, "", 1, features[i])]
  }

  printf "  " features[2] " -> "
  if (n == 3)  {
    printf features[3] " [color=\"" color[2] "\"] " 
  } else {
    printf "{"
    for (i=3; i<=n; i++) {
      printf features[i] " "
    }
    printf "}"
  }

  print " " depstyle[deptype] ";"
}

