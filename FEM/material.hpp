#pragma once

#include <iostream>
#include <string>
#include <vector>

class Material
{
public:
	Material();
	Material(std::string name);

	class Properties
	{
	public:
		double GetDensity() const
		{
			return density;
		}
		void SetDensity(double newDensity)
		{
			density = newDensity;
		}

		friend std::istream& operator>>(std::istream& is, Properties& prop);

	private:
		double density;
	};

	const std::string& GetName() const
	{
		return name;
	}

	const Properties& GetProperties() const
	{
		return properties;
	}
	void SetProperties(const Properties& prop)
	{
		properties = prop;
	}

	friend std::istream& operator>>(std::istream& is, Material& material);

private:
	std::string name;
	Properties properties;
};

inline bool operator==(const Material& lhs, const Material& rhs)
{
	return lhs.GetName() == rhs.GetName();
}

namespace std
{
template <>
struct hash<Material>
{
	std::size_t operator()(Material const& material) const noexcept
	{
		std::size_t h = std::hash<std::string>{}(material.GetName());
		return h;
	}
};
}  // namespace std
