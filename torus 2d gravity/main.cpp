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

#define numGravityFields 2 // max is torusPowerRange - 2, anything greater is experimental, min 1, greater then torusPowerRange should give no better results
#define gravityConstant 1

#define torusPowerRange 7 // max is 32, min is 3

#define halfScreenx 600
#define halfScreeny 400

#define massFactor 3.14
#define minRadius 0.05 // can be anything less then maxRadius
#define maxRadius 0.5 // max is 0.5 because diamiter is 1 and space is partitioned every unit

#define numBalls 10000

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
	uint32_t m_z;
	uint32_t m_w;

	uint32_t torusRange;
	uint32_t torusMod;

	vector<ball*> balls;

	unordered_map<uint32_t, unordered_map<uint32_t, vector<ball*>>> collisionSpace;

	unordered_map<uint32_t, unordered_map<uint32_t, vector<ball*>>> renderSpace;

	unordered_map<uint32_t, unordered_map<uint32_t, clusterBall>> layeredGravityFields[numGravityFields];

	uint32_t intRand()
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

	void ballPullBall(clusterBall* ball1, clusterBall* ball2, int x, int y, uint32_t unitDistance)
	{
		vd2d dpos = vd2d{ double(x), double(y) } + ball2->pos - ball2->pos.touint() + ball1->pos.touint() - ball1->pos;
		double dis = dpos.mag2();

		dpos *= gravityConstant / (dis * sqrt(dis));

		ball1->posv += dpos * ball2->mass;
		ball2->posv -= dpos * ball1->mass;
	}

	void gravity()
	{
		for (int i = 0; i < balls.size(); i++)
		{
			layeredGravityFields[0][uint32_t(balls[i]->pos.x)][uint32_t(balls[i]->pos.y)].mass += balls[i]->mass;
			layeredGravityFields[0][uint32_t(balls[i]->pos.x)][uint32_t(balls[i]->pos.y)].pos += balls[i]->pos * balls[i]->mass;
		}

		for (int i = 0; i < numGravityFields - 1; i++)
		{
			for (auto j = layeredGravityFields[i].begin(); j != layeredGravityFields[i].end(); j++)
			{
				for (auto k = j->second.begin(); k != j->second.end(); k++)
				{
					layeredGravityFields[i + 1][j->first >> 1][k->first >> 1].mass += k->second.mass;
					layeredGravityFields[i + 1][j->first >> 1][k->first >> 1].pos += k->second.pos;
					k->second.pos /= k->second.mass;
				}
			}
		}

		for (auto i = layeredGravityFields[numGravityFields - 1].begin(); i != layeredGravityFields[numGravityFields - 1].end(); i++)
		{
			for (auto j = i->second.begin(); j != i->second.end(); j++)
			{
				j->second.pos /= j->second.mass;
			}
		}

		unordered_map<uint32_t, unordered_map<uint32_t, clusterBall>>::iterator findx;
		unordered_map<uint32_t, clusterBall>::iterator findy;

		for (int i = numGravityFields - 1; i >= 0; i--)
		{
			uint32_t unitDistance = 1 << i;

			for (auto j = layeredGravityFields[i].begin(); j != layeredGravityFields[i].end(); j++)
			{
				for (auto k = j->second.begin(); k != j->second.end(); k++)
				{
					uint32_t anchorx, anchory;
					int unitDifx, unitDify;
					anchorx = j->first & 0xfffffffe;
					anchory = k->first & 0xfffffffe;

					for (int x = 0; x < 4; x++)
					{
						findx = layeredGravityFields[i].find(uint32_t(j->first + x) & (torusMod >> i));

						if (findx != layeredGravityFields[i].end())
						{
							for (int y = -2; y < 4; y++)
							{
								if (x > 1 || y > 1)
								{
									unitDifx = anchorx + x - j->first;
									unitDify = anchory + y - k->first;
									if (abs(unitDifx) > 1 || abs(unitDify) > 1)
									{
										findy = findx->second.find(uint32_t(k->first + y) & (torusMod >> i));

										if (findy != findx->second.end())
										{
											ballPullBall(&k->second, &findy->second, unitDifx, unitDify, unitDistance);
										}
									}
								}
							}
						}
					}

					if (i != 0)
					{
						for (int x = 0; x < 2; x++)
						{
							findx = layeredGravityFields[i - 1].find(j->first << 1 & x);

							if (findx != layeredGravityFields[i - 1].end())
							{
								for (int y = 0; y < 2; y++)
								{
									findy = findx->second.find(k->first << 1 & y);

									if (findy != findx->second.end())
									{
										findy->second.posv += k->second.pos;
									}
								}
							}
						}
					}
				}
			}
		}

		for (int i = 0; i < balls.size(); i++)
		{
			balls[i]->posv += layeredGravityFields[0][uint32_t(balls[i]->pos.x)][uint32_t(balls[i]->pos.y)].posv;
		}

		for (int i = 0; i < numGravityFields; i++)
		{
			layeredGravityFields[i].clear();
		}
	}

	void ballToBall(ball* ball1, ball* ball2, int x, int y)
	{
		vd2d dpos = vd2d{ double(x), double(y) } + ball2->pos - ball2->pos.touint() + ball1->pos.touint() - ball1->pos;
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
		for (int i = 0; i < balls.size(); i++)
			collisionSpace[uint32_t(balls[i]->pos.x)][uint32_t(balls[i]->pos.y)].push_back(balls[i]);

		unordered_map<uint32_t, unordered_map<uint32_t, vector<ball*>>>::iterator findx;
		unordered_map<uint32_t, vector<ball*>>::iterator findy;

		for (auto partitionx = collisionSpace.begin(); partitionx != collisionSpace.end(); partitionx++)
		{
			for (auto partitiony = partitionx->second.begin(); partitiony != partitionx->second.end(); partitiony++)
			{
				for (int b1 = 0; b1 < partitiony->second.size(); b1++)
				{
					for (int b2 = b1 + 1; b2 < partitiony->second.size(); b2++)
						ballToBall(partitiony->second[b1], partitiony->second[b2], 0, 0);

					findy = partitionx->second.find(uint32_t(partitiony->first + 1) & torusMod);

					if (findy != partitionx->second.end())
						for (int b2 = 0; b2 < findy->second.size(); b2++)
							ballToBall(partitiony->second[b1], findy->second[b2], 0, 1);

					findx = collisionSpace.find(uint32_t(partitionx->first + 1) & torusMod);

					if (findx != collisionSpace.end())
					{
						findy = findx->second.find(uint32_t(partitiony->first + 1) & torusMod);

						if (findy != findx->second.end())
							for (int b2 = 0; b2 < findy->second.size(); b2++)
								ballToBall(partitiony->second[b1], findy->second[b2], 1, 1);

						findy = findx->second.find(partitiony->first);

						if (findy != findx->second.end())
							for (int b2 = 0; b2 < findy->second.size(); b2++)
								ballToBall(partitiony->second[b1], findy->second[b2], 1, 0);

						findy = findx->second.find(uint32_t(partitiony->first - 1) & torusMod);

						if (findy != findx->second.end())
							for (int b2 = 0; b2 < findy->second.size(); b2++)
								ballToBall(partitiony->second[b1], findy->second[b2], 1, -1);
					}
				}
			}
		}

		collisionSpace.clear();
	}

	void moveBalls(float fElapsedTime)
	{
		for (int i = 0; i < balls.size(); i++)
		{
			balls[i]->pos += balls[i]->posv * fElapsedTime;
			balls[i]->pos = balls[i]->pos - (balls[i]->pos / torusRange).floor() * torusRange;
		}
	}

	void drawBalls()
	{
		Clear(Pixel(0, 0, 0));

		for (int i = 0; i < balls.size(); i++)
			renderSpace[uint32_t(balls[i]->pos.x)][uint32_t(balls[i]->pos.y)].push_back(balls[i]);

		vd2d startPos = (pos - halfScreen / zoom).floor();
		vd2d endPos = (pos + halfScreen / zoom).floor();

		unordered_map<uint32_t, unordered_map<uint32_t, vector<ball*>>>::iterator findx;
		unordered_map<uint32_t, vector<ball*>>::iterator findy;

		for (int partitionx = startPos.x; partitionx <= endPos.x; partitionx++)
		{
			findx = renderSpace.find(uint32_t(partitionx) & torusMod);

			if (findx != renderSpace.end())
			{
				for (int partitiony = startPos.y; partitiony <= endPos.y; partitiony++)
				{
					findy = findx->second.find(uint32_t(partitiony) & torusMod);

					if (findy != findx->second.end())
					{
						for (int b = 0; b < findy->second.size(); b++)
						{
							vd2d bPos = vd2d{ double(partitionx), double(partitiony) } + findy->second[b]->pos - findy->second[b]->pos.touint();

							FillCircle((bPos - pos) * zoom + halfScreen, zoom * findy->second[b]->radius, findy->second[b]->color);
						}
					}
				}
			}
		}

		renderSpace.clear();
	}

	bool OnUserCreate() override
	{
		pos = { 0,0 };
		zoom = 64;

		halfScreen = { halfScreenx, halfScreeny };
		m_z = (uint32_t)duration_cast<seconds>(high_resolution_clock::now().time_since_epoch()).count();
		m_w = (uint32_t)duration_cast<microseconds>(high_resolution_clock::now().time_since_epoch()).count();

		torusRange = 1 << torusPowerRange;
		torusMod = 0xffffffff >> (32 - torusPowerRange);

		for (int i = 0; i < numBalls; i++)
		{
			double randNum = doubleRand() * 6.28318530718;
			vd2d bPos = (vd2d{ doubleRand(), doubleRand() }) * torusRange;
			vd2d bPosv = /*vd2d{ 0, 0 };//*/ (vd2d{ cos(randNum), sin(randNum) }) * 2;
			Pixel bColor = mapToRainbow(doubleRand() * 6);
			balls.push_back(new ball(bPos, bPosv, bColor, doubleRand() * (maxRadius - minRadius) + minRadius));
		}

		return true;
	}

	bool OnUserUpdate(double fElapsedTime) override
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