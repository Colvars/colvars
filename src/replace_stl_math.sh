for f in *.{cpp,h} ; do

    # The following is not the full list under std::, but should catch all
    for fn in cos sin tan acos asin atan atan2 exp erf erfc log10 floor ceil fabs round sqrt ; do
        sed -i "s/std::${fn}/cvm::${fn}/g" $f
    done

    # Special cases
    sed -i 's/std::abs/cvm::fabs/g' $f
    sed -i 's/std::log/cvm::logn/g' $f

    # NOTE: Some calls should be changed to cvm::power_integer()
    sed -i 's/std::pow/cvm::pow/g' $f

done
