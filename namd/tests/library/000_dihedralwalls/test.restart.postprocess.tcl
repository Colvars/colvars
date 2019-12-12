exec sed -n "/restraint/q;p" test.colvars.state > /tmp/test.colvars.state
exec mv -f /tmp/test.colvars.state test.colvars.state

exec sed -n "/restraint/q;p" test.restart.colvars.state > /tmp/test.restart.colvars.state
exec mv -f /tmp/test.restart.colvars.state test.restart.colvars.state
