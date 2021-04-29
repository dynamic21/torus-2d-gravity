#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"
#include <unordered_map>

using std::max;
using std::min;
using std::endl;
using std::cout;
using std::vector;
using std::unordered_map;

using olc::Key;
using olc::vd2d;
using olc::Pixel;

using std::chrono::seconds;
using std::chrono::microseconds;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;

#define gravityConstant 1
#define massFactor 3.14
#define minRadius 0.05
#define maxRadius 0.5
#define torusRange 16 // 2 ^ 8
#define torusPowerRange 4 // 2 ^ 8
#define halfScreenx 600
#define halfScreeny 400

class ball
{
public:
	vd2d pos;
	vd2d posv;
	Pixel color;
	double mass;
	double radius;

	ball(vd2d Pos, vd2d Posv, Pixel Color, double Radius)
	{
		pos = Pos;
		posv = Posv;
		color = Color;
		radius = Radius;
		mass = massFactor * radius * radius;
	}

	/*ball(ball* Ball)
	{
		pos = Ball->pos - vd2d{ floor(Ball->pos.x), floor(Ball->pos.y) };
		posv = Ball->posv;
		color = Ball->color;
		radius = Ball->radius;
		mass = Ball->mass;
	}*/
};

class clusterBall
{
public:
	vd2d pos;
	vd2d posv;
	double mass;

	clusterBall(vd2d Pos, vd2d Posv, double Mass)
	{
		pos = Pos;
		posv = Posv;
		mass = Mass;
	}
};

class Example : public olc::PixelGameEngine
{
public:
	Example()
	{
		sAppName = "2D Torus Universe";
	}

public:
	vd2d pos;
	double zoom;

	vd2d halfScreen;
	unsigned int m_z;
	unsigned int m_w;

	vector<ball*> balls;

	unordered_map<unsigned int, unordered_map<unsigned int, vector<ball*>>> grid;

	unsigned int intRand()
	{
		m_z = 36969 * (m_z & 65535) + (m_z >> 16);
		m_w = 18000 * (m_w & 65535) + (m_w >> 16);

		return (m_z << 16) + m_w;
	}

	double doubleRand() { return (intRand() + 1.0) * 2.328306435454494e-10; }

	Pixel mapToRainbow(double d)
	{
		double r = (d > 3) ? max(0.0, min(1.0, d - 4)) : max(0.0, min(1.0, 2 - d));
		double g = (d > 2) ? max(0.0, min(1.0, 4 - d)) : max(0.0, min(1.0, d));
		double b = (d > 4) ? max(0.0, min(1.0, 6 - d)) : max(0.0, min(1.0, d - 2));

		return Pixel(r * 0xff, g * 0xff, b * 0xff);
	}

	void userControl(double fElapsedTime)
	{
		if (GetKey(Key::Q).bHeld) { zoom /= pow(2, fElapsedTime); }
		if (GetKey(Key::E).bHeld) { zoom *= pow(2, fElapsedTime); }

		if (GetKey(Key::W).bHeld || GetKey(Key::UP).bHeld) { pos.y -= 200 * fElapsedTime / zoom; }
		if (GetKey(Key::A).bHeld || GetKey(Key::LEFT).bHeld) { pos.x -= 200 * fElapsedTime / zoom; }
		if (GetKey(Key::S).bHeld || GetKey(Key::DOWN).bHeld) { pos.y += 200 * fElapsedTime / zoom; }
		if (GetKey(Key::D).bHeld || GetKey(Key::RIGHT).bHeld) { pos.x += 200 * fElapsedTime / zoom; }
	}

	/*void ballToBall(int i, int j)
	{
		vd2d dpos = balls[j].pos - balls[i].pos;
		double dis = dpos.mag2();

		if (dis < 1)
		{
			dpos /= sqrt(dis);
			dis = (balls[j].posv - balls[i].posv).dot(dpos);

			if (dis < 0)
			{
				dpos *= dis;
				balls[i].posv += dpos;
				balls[j].posv -= dpos;
			}
		}
	}*/

	void collision()
	{
		unordered_map<unsigned int, unordered_map<unsigned int, vector<ball*>>>::iterator findx;
		unordered_map<unsigned int, vector<ball*>>::iterator findy;

		for (auto i = grid.begin(); i != grid.end(); i++)
		{
			for (auto j = i->second.begin(); j != i->second.end(); j++)
			{
				for (int k = 0; k < j->second.size(); k++)
				{
					for (int x = 0; x < 2; x++)
					{
						findx = grid.find(i->first + x);

						if (findx != grid.end())
						{
							for (int y = -1; y < 2; y++)
							{
								if (x != 0 && y != -1)
								{
									findy = findx->second.find(j->first + y);

									if (findy != findx->second.end())
									{
										for (int l = 0; l < findy->second.size(); l++)
										{
											//ballToBall(j->second[k], findy->second[l]);
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	void moveBalls(float fElapsedTime)
	{
		for (int i = 0; i < balls.size(); i++)
		{
			balls[i]->pos += balls[i]->posv * fElapsedTime;
			balls[i]->pos = balls[i]->pos - torusRange * (balls[i]->pos / torusRange).floor();
		}
	}

	void drawBalls()
	{
		Clear(Pixel(0, 0, 0));

		for (int i = 0; i < balls.size(); i++)
		{
			grid[floor(balls[i]->pos.x)][floor(balls[i]->pos.y)].push_back(balls[i]);
		}

		/*vd2d roundPos = vd2d{ unsigned int(pos.x) - double(pos.x < 0), unsigned int(pos.y) - double(pos.y < 0) };
		vd2d modPos = roundPos - pos;
		vd2d cornerPos = halfScreen / zoom + modPos;
		cornerPos = vd2d{ unsigned int(cornerPos.x) - double(cornerPos.x < 0), unsigned int(cornerPos.y) - double(cornerPos.y < 0) };

		for (int i = 0; i < balls.size(); i++)
		{
			FillCircle((balls[i]->pos - pos) * zoom + halfScreen, zoom * balls[i]->radius, balls[i]->color);
		}*/
		vd2d startPos = (pos - halfScreen / zoom).floor();
		vd2d endPos = (pos + halfScreen / zoom).floor();

		unordered_map<unsigned int, unordered_map<unsigned int, vector<ball*>>>::iterator findx;
		unordered_map<unsigned int, vector<ball*>>::iterator findy;

		for (int mx = startPos.x; mx < endPos.x; mx++)
		{
			findx = grid.find(mx);

			if (findx != grid.end())
			{
				for (int my = startPos.y; my < endPos.y; my++)
				{
					findy = findx->second.find(my);

					if (findy != findx->second.end())
					{
						for (int i = 0; i < grid[mx][my].size(); i++)
						{
							vd2d bPos = vd2d{ double(mx), double(my) } + grid[mx][my][i]->pos - grid[mx][my][i]->pos.floor();

							FillCircle((bPos - pos) * zoom + halfScreen, zoom * grid[mx][my][i]->radius, grid[mx][my][i]->color);
						}
					}
				}
			}
		}
	}

	bool OnUserCreate() override
	{
		pos = { 0,0 };
		zoom = 16;

		halfScreen = { halfScreenx, halfScreeny };
		m_z = (unsigned int)duration_cast<seconds>(high_resolution_clock::now().time_since_epoch()).count();
		m_w = (unsigned int)duration_cast<microseconds>(high_resolution_clock::now().time_since_epoch()).count();

		for (int i = 0; i < 10; i++)
		{
			double randNum = doubleRand() * 6.28318530718;
			vd2d bPos = (vd2d{ doubleRand(), doubleRand() }) * pow(2, torusPowerRange);
			vd2d bPosv = (vd2d{ cos(randNum), sin(randNum) }) * 1;
			Pixel bColor = mapToRainbow(doubleRand() * 6);
			balls.push_back(new ball(bPos, bPosv, bColor, doubleRand() * (maxRadius - minRadius) + minRadius));
		}

		return true;
	}

	bool OnUserUpdate(float fElapsedTime) override
	{
		userControl(fElapsedTime);
		//gravity();
		//collision();
		moveBalls(fElapsedTime);
		drawBalls();

		return true;
	}
};

int main()
{
	Example demo;

	if (demo.Construct(halfScreenx * 2, halfScreeny * 2, 1, 1))
		demo.Start();

	return 0;
}