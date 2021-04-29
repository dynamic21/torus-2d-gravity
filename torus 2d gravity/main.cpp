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
#define minRadius 0.05 // can be anything less then 0.5
#define maxRadius 0.5 // must be 0.5 cuz diamiter is 1 and collision partitions by 1
#define torusRange 8 // 2 ^ 3
#define torusPowerRange 3 // 2 ^ 3
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
};

class clusterBall
{
public:
	vd2d pos;
	vd2d posv;
	double mass;

	clusterBall()
	{
		pos = { 0,0 };
		posv = { 0,0 };
		mass = 0;
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

	unordered_map<unsigned int, unordered_map<unsigned int, vector<ball*>>> partitionedSpace;

	unordered_map<unsigned int, unordered_map<unsigned int, clusterBall*>> layeredGravityFields;

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

		if (GetKey(Key::W).bHeld || GetKey(Key::UP).bHeld) { pos.y -= 400 * fElapsedTime / zoom; }
		if (GetKey(Key::A).bHeld || GetKey(Key::LEFT).bHeld) { pos.x -= 400 * fElapsedTime / zoom; }
		if (GetKey(Key::S).bHeld || GetKey(Key::DOWN).bHeld) { pos.y += 400 * fElapsedTime / zoom; }
		if (GetKey(Key::D).bHeld || GetKey(Key::RIGHT).bHeld) { pos.x += 400 * fElapsedTime / zoom; }
	}

	void ballPullBall(clusterBall* ball1, clusterBall* ball2, int mx, int my)//add unit distance
	{
		vd2d dpos = vd2d{ double(mx), double(my) } + ball2->pos - ball2->pos.floor() + ball1->pos.floor() - ball1->pos;
		double dis = dpos.mag2();

		dpos *= gravityConstant / (dis * sqrt(dis));

		ball1->posv += dpos * ball2->mass;
		ball2->posv -= dpos * ball1->mass;
	}

	void gravity()
	{
		//
	}

	void ballToBall(ball* ball1, ball* ball2, int mx, int my)
	{
		vd2d dpos = vd2d{ double(mx), double(my) } + ball2->pos - ball2->pos.floor() + ball1->pos.floor() - ball1->pos;
		double dis = dpos.mag();

		if (dis < ball2->radius + ball1->radius)
		{
			dpos /= dis;
			dis = (ball2->posv - ball1->posv).dot(dpos);

			if (dis < 0)
			{
				dpos *= 2 * dis / (ball1->mass + ball2->mass);
				ball1->posv += dpos * ball2->mass;
				ball2->posv -= dpos * ball1->mass;
			}
		}
	}

	void collision()
	{
		unordered_map<unsigned int, unordered_map<unsigned int, vector<ball*>>>::iterator findx;
		unordered_map<unsigned int, vector<ball*>>::iterator findy;

		unsigned int tx, ty;

		for (auto i = partitionedSpace.begin(); i != partitionedSpace.end(); i++)
		{
			for (auto j = i->second.begin(); j != i->second.end(); j++)
			{
				for (int k = 0; k < j->second.size(); k++)
				{
					for (int l = k + 1; l < j->second.size(); l++)
					{
						ballToBall(j->second[k], j->second[l], 0, 0);
					}

					for (int mx = 0; mx < 2; mx++)
					{
						tx = unsigned int(i->first + mx) % torusRange;
						findx = partitionedSpace.find(tx);

						if (findx != partitionedSpace.end())
						{
							for (int my = -1; my < 2; my++)
							{
								if (mx > 0 || my > 0)
								{
									ty = unsigned int(j->first + my) % torusRange;
									findy = findx->second.find(ty);

									if (findy != findx->second.end())
									{
										for (int l = 0; l < findy->second.size(); l++)
										{
											ballToBall(j->second[k], findy->second[l], mx, my);
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

		partitionedSpace.clear();

		for (int i = 0; i < balls.size(); i++)
		{
			partitionedSpace[floor(balls[i]->pos.x)][floor(balls[i]->pos.y)].push_back(balls[i]);
		}

		vd2d startPos = (pos - halfScreen / zoom).floor();
		vd2d endPos = vd2d{ 1.0, 1.0 } + (pos + halfScreen / zoom).floor();

		unordered_map<unsigned int, unordered_map<unsigned int, vector<ball*>>>::iterator findx;
		unordered_map<unsigned int, vector<ball*>>::iterator findy;

		unsigned int tx, ty;

		for (int mx = startPos.x; mx < endPos.x; mx++)
		{
			tx = unsigned int(mx) % torusRange;
			findx = partitionedSpace.find(tx);

			if (findx != partitionedSpace.end())
			{
				for (int my = startPos.y; my < endPos.y; my++)
				{
					ty = unsigned int(my) % torusRange;
					findy = findx->second.find(ty);

					if (findy != findx->second.end())
					{
						for (int i = 0; i < partitionedSpace[tx][ty].size(); i++)
						{
							vd2d bPos = vd2d{ double(mx), double(my) } + partitionedSpace[tx][ty][i]->pos - partitionedSpace[tx][ty][i]->pos.floor();

							FillCircle((bPos - pos) * zoom + halfScreen, zoom * partitionedSpace[tx][ty][i]->radius, partitionedSpace[tx][ty][i]->color);
						}
					}
				}
			}
		}
	}

	bool OnUserCreate() override
	{
		pos = { 0,0 };
		zoom = 64;

		halfScreen = { halfScreenx, halfScreeny };
		m_z = (unsigned int)duration_cast<seconds>(high_resolution_clock::now().time_since_epoch()).count();
		m_w = (unsigned int)duration_cast<microseconds>(high_resolution_clock::now().time_since_epoch()).count();

		for (int i = 0; i < 10; i++)
		{
			double randNum = doubleRand() * 6.28318530718;
			vd2d bPos = (vd2d{ doubleRand(), doubleRand() }) * torusRange;
			vd2d bPosv = (vd2d{ cos(randNum), sin(randNum) }) * 2;
			Pixel bColor = mapToRainbow(doubleRand() * 6);
			balls.push_back(new ball(bPos, bPosv, bColor, doubleRand() * (maxRadius - minRadius) + minRadius));
		}

		return true;
	}

	bool OnUserUpdate(float fElapsedTime) override
	{
		userControl(fElapsedTime);
		//gravity();
		collision();
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