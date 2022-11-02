#pragma once
#include <stddef.h>

namespace acamcad {
/**
* 两个颜色和纹理坐标接近小于这个值，会被认为是同一个点
*/
const float ColorThreshold = 0.001f;

/**
* The RGB is 0-255, 不允许输入小于0的数
* The RGBf is 0-1; 小于0的数会被设置为0；
*/

struct MRGB
{
	uint32_t r, g, b;
	MRGB(uint32_t red, uint32_t green, uint32_t blue)
	{
		r = red % 255;
		g = green % 255;
		b = blue % 255;
	}
	MRGB() : r(0), g(0), b(0) {}

	uint32_t setR(uint32_t red) { return r = red % 255; }
	uint32_t setG(uint32_t green) { return g = green % 255; }
	uint32_t setB(uint32_t blue) { return b = blue % 255; }

	uint32_t R() { return r; }
	uint32_t G() { return g; }
	uint32_t B() { return b; }
};

struct MRGBf
{
	float r, g, b;
	MRGBf(float red, float green, float blue)
	{
		red = red >= 0 ? red : 0;
		r = red - (int)red;
		green = green >= 0 ? green : 0;
		g = green - (int)green;
		blue = blue >= 0 ? blue : 0;
		b = blue - (int)blue;
	}
	MRGBf() : r(0), g(0), b(0) {}

	float setR(float red) {
		red = red >= 0 ? red : 0;
		return r = red - (int)red;
	}
	float setG(float green) {
		green = green >= 0 ? green : 0;
		return g = green - (int)green;
	}
	float setB(float blue) {
		blue = blue >= 0 ? blue : 0;
		return b = blue - (int)blue;
	}

	float R() { return r; }
	float G() { return g; }
	float B() { return b; }
};

struct MRGBA
{
	uint32_t r, g, b, a;

	MRGBA(uint32_t red, uint32_t green, uint32_t blue, uint32_t alpha)
	{
		r = red % 255;
		g = green % 255;
		b = blue % 255;

		if (alpha > 100) alpha = 100;
		a = alpha;
	}

	MRGBA(uint32_t red, uint32_t green, uint32_t blue) : a(0)
	{
		r = red % 255;
		g = green % 255;
		b = blue % 255;
	}

	uint32_t setR(uint32_t red) { return r = red % 255; }
	uint32_t setG(uint32_t green) { return g = green % 255; }
	uint32_t setB(uint32_t blue) { return b = blue % 255; }
	void setA(uint32_t alpha) {
		if (alpha > 100) alpha = 100;
		a = alpha;
	}

	uint32_t R() { return r; }
	uint32_t G() { return g; }
	uint32_t B() { return b; }
	float A() { return a / 100.0f; }
};

struct MRGBAf
{
	float r, g, b, a;
	MRGBAf(float red, float green, float blue, float alpha)
	{
		red = red >= 0 ? red : 0;
		r = red - (int)red;
		green = green >= 0 ? green : 0;
		g = green - (int)green;
		blue = blue >= 0 ? blue : 0;
		b = blue - (int)blue;

		if (alpha < 0) alpha = 0;
		if (alpha > 1) alpha = 1;
		a = alpha;
	}

	MRGBAf() : r(0), g(0), b(0), a(0) {}

	float setR(float red) {
		red = red >= 0 ? red : 0;
		return r = red - (int)red;
	}
	float setG(float green) {
		green = green >= 0 ? green : 0;
		return g = green - (int)green;
	}
	float setB(float blue) {
		blue = blue >= 0 ? blue : 0;
		return b = blue - (int)blue;
	}
	float setA(float alpha) {
		if (alpha < 0) alpha = 0;
		if (alpha > 1) alpha = 1;
		return a = alpha;
	}

	float R() { return r; }
	float G() { return g; }
	float B() { return b; }
	float A() { return a; }
};

struct Texcoord
{
	float uv[3];

	Texcoord()
	{
		uv[0] = 0; uv[1] = 1; uv[2] = 0;
	}

	Texcoord(float u, float v)
	{
		uv[0] = u; uv[1] = v; uv[2] = 0;
	}
	Texcoord(float u, float v, float w)
	{
		uv[0] = u; uv[1] = v; uv[2] = w;
	}

	void setW(float w) { uv[2] = w; }

	float& operator[](unsigned i) { return uv[i]; } 
	float operator[](unsigned i) const { return uv[i]; }

	bool operator<(const Texcoord &rhs) const {
		
		for (int i = 0; i < 3; i++)
		{
			if (uv[i] < rhs.uv[i] - ColorThreshold)
				return true;
			else if (uv[i] > rhs.uv[i] + ColorThreshold)
				return false;
		}
		return true;
	}
};

inline MRGBf RGB_to_RGBF(MRGB& color)
{
	return(MRGBf(color.r / 255.0f, color.g / 255.0f, color.b / 255.0f));
}

inline MRGB RGBF_to_RGB(MRGBf& color)
{
	return(MRGB(uint32_t(color.r * 255), uint32_t(color.g * 255), uint32_t(color.b * 255)));
}

inline MRGBAf RGBA_to_RGBAF(MRGBA& color)
{
	return(MRGBAf(color.r / 255.0f, color.g / 255.0f, color.b / 255.0f, color.a / 100.0f));
}

inline MRGBA RGBAF_to_RGBA(MRGBAf& color)
{
	return(MRGBA(uint32_t(color.r * 255), uint32_t(color.g * 255), uint32_t(color.b * 255), uint32_t(color.a * 100)));
}

/**
* Determine if the colors are equal
*/
inline bool equality(Texcoord& lhs, Texcoord& rhs)
{
	bool t0 = (lhs[0] - rhs[0] < ColorThreshold) && (rhs[0] - lhs[0] < ColorThreshold);
	bool t1 = (lhs[1] - rhs[1] < ColorThreshold) && (rhs[1] - lhs[1] < ColorThreshold);
	bool t2 = (lhs[2] - rhs[2] < ColorThreshold) && (rhs[2] - lhs[2] < ColorThreshold);

	if (t0 && t1 && t2)
	{
		return true;
	}
	else
	{
		return false;
	}
}

inline bool equality(MRGB& lhs, MRGB& rhs)
{
	if ((lhs.r - rhs.r) == 0 && (lhs.g - rhs.g) == 0 && (lhs.b - rhs.b) == 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

inline bool equality(MRGBf& lhs, MRGBf& rhs)
{
	bool r = (lhs.r - rhs.r < ColorThreshold) && (rhs.r - lhs.r < ColorThreshold);
	bool g = (lhs.g - rhs.g < ColorThreshold) && (rhs.g - lhs.g < ColorThreshold);
	bool b = (lhs.b - rhs.b < ColorThreshold) && (rhs.b - lhs.b < ColorThreshold);

	if (r && g && b)
	{
		return true;
	}
	else
	{
		return false;
	}
}

inline bool equality(MRGBAf& lhs, MRGBAf& rhs)
{
	bool r = (lhs.r - rhs.r < ColorThreshold) && (rhs.r - lhs.r < ColorThreshold);
	bool g = (lhs.g - rhs.g < ColorThreshold) && (rhs.g - lhs.g < ColorThreshold);
	bool b = (lhs.b - rhs.b < ColorThreshold) && (rhs.b - lhs.b < ColorThreshold);
	bool a = (lhs.a - rhs.a < ColorThreshold) && (rhs.a - lhs.a < ColorThreshold);

	if (r && g && b && a)
	{
		return true;
	}
	else
	{
		return false;
	}
}


} //namespace acamcad
