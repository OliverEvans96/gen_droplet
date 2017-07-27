#!/bin/bash

# Generate water droplet, shifted to zero

radius=70
shape="sphere"

echo ".x waterdroplet_tip4p_new.C(${radius}, \"${shape}\");" | root -lbq

echo
echo "Generated 'waterdroplet.shifted.to.zero.spce.${radius}'"
