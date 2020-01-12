#include "stdafx.hpp"

#include "material.hpp"

std::istream& operator>>(std::istream& is, Material::Properties& prop)
{
	is >> prop.density;
	return is;
}

/*
 * example of input:
 * Al_Alloy 0.0028
 */
std::istream& operator>>(std::istream& is, Material& material)
{
	is >> material.name >> material.properties;
	return is;
}

Material::Material()
{
}

Material::Material(std::string name) : name(std::move(name))
{
}
