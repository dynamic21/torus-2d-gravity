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
#define gravityConstant 10

#define torusPowerRange 4 // max is 32, min is 3

#define halfScreenx 600
#define halfScreeny 400

#define massFactor 3.14
#define minRadius 0.05 // can be anything less then maxRadius
#define maxRadius 0.5 // max is 0.5 because diamiter is 1 and space is partitioned every unit

#define numBalls 2

class ball
{
public:
	vd2d pos;
	vd2d vel;
	Pixel color;
	double mass;
	double radius;

	ball(vd2d Pos, vd2d Posv, Pixel Color, double Radius)
	{
		pos = Pos;
		vel = Posv;
		color = Color;
		radius = Radius;
		mass = massFactor * radius * radius;
	}
};

class clusterBall
{
public:
	vd2d pos;
	vd2d vel;
	double mass;

	clusterBall()
	{
		pos = vd2d(0, 0);
		vel = vd2d(0, 0);
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

	vector<ball*> balls = { new ball(vd2d(10, 30), vd2d(0, 0), Pixel(0xff, 0, 0xff), 0.5), new ball(vd2d(10, 10), vd2d(0, 0), Pixel(0xff, 0xff, 0), 0.5) };

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

	void ballPullBall(clusterBall* ball1, clusterBall* ball2, int dx, int dy, uint32_t unitDistance, uint32_t partion1x, uint32_t partion1y, uint32_t partion2x, uint32_t partion2y)
	{
		vd2d dpos = vd2d(double(dx) - partion2x + partion1x, double(dy) - partion2y + partion1y) * unitDistance + ball2->pos - ball1->pos;
		double dis = dpos.mag2();

		dpos *= gravityConstant / (dis * sqrt(dis));

		ball1->vel += dpos * ball2->mass;
		//cout << "before: " << ball2->vel << ", after:" << ball2->vel - dpos * ball1->mass << endl;
		ball2->vel -= dpos * ball1->mass;
	}

	void gravity(double fElapsedTime)
	{
		for (int b = 0; b < balls.size(); b++)
		{
			layeredGravityFields[0][uint32_t(balls[b]->pos.x)][uint32_t(balls[b]->pos.y)].mass += balls[b]->mass;
			layeredGravityFields[0][uint32_t(balls[b]->pos.x)][uint32_t(balls[b]->pos.y)].pos += balls[b]->pos * balls[b]->mass;
		}

		for (int layer = 0; layer < numGravityFields - 1; layer++)
		{
			for (auto fieldx = layeredGravityFields[layer].begin(); fieldx != layeredGravityFields[layer].end(); fieldx++)
			{
				for (auto fieldy = fieldx->second.begin(); fieldy != fieldx->second.end(); fieldy++)
				{
					uint32_t x = fieldx->first >> 1;
					uint32_t y = fieldy->first >> 1;

					layeredGravityFields[layer + 1][x][y].mass += fieldy->second.mass;
					layeredGravityFields[layer + 1][x][y].pos += fieldy->second.pos;
					fieldy->second.pos /= fieldy->second.mass;
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

		for (int layer = numGravityFields - 1; layer >= 0; layer--)
		{
			uint32_t unitDistance = 1 << layer;
			uint32_t layeredTorusMod = torusMod >> layer;

			for (auto fieldx = layeredGravityFields[layer].begin(); fieldx != layeredGravityFields[layer].end(); fieldx++)
			{
				for (auto fieldy = fieldx->second.begin(); fieldy != fieldx->second.end(); fieldy++)
				{
					int unitDifx, unitDify;

					uint32_t anchorx = fieldx->first & 0xfffffffe;
					uint32_t anchory = fieldy->first & 0xfffffffe;

					/*for (int x = 0; x < 4; x++)
					{
						int tx = anchorx + x & (torusMod >> layer);

						for (int y = -2; y < 4; y++)
						{
							int ty = anchory + y & (torusMod >> layer);

							if (x > 1 || y > 1)
							{
								unitDifx = anchorx - fieldx->first + x;
								unitDify = anchory - fieldy->first + y;

								if (abs(unitDifx) > 1 || abs(unitDify) > 1)
								{
									vd2d bPos = vd2d(double(tx + 0.5), double(ty + 0.5)) * unitDistance;
									FillCircle((bPos - pos) * zoom + halfScreen, zoom * unitDistance / 2, Pixel(0, 0, (numGravityFields - layer) * (250 / numGravityFields)));
								}
							}
						}
					}*/

					for (int dx = 0; dx < 4; dx++)
					{
						findx = layeredGravityFields[layer].find(anchorx + dx & layeredTorusMod);

						if (findx != layeredGravityFields[layer].end())
						{
							for (int dy = -2; dy < 4; dy++)
							{
								if (dx > 1 || dy > 1)
								{
									unitDifx = anchorx - fieldx->first + dx;
									unitDify = anchory - fieldy->first + dy;

									if (abs(unitDifx) > 1 || abs(unitDify) > 1)
									{
										findy = findx->second.find(anchory + dy & layeredTorusMod);

										if (findy != findx->second.end())
										{
											ballPullBall(&fieldy->second, &findy->second, unitDifx, unitDify, unitDistance, fieldx->first, fieldy->first, findx->first, findy->first);
										}
									}
								}
							}
						}
					}

					for (auto i = layeredGravityFields[1].begin(); i != layeredGravityFields[1].end(); i++)
					{
						for (auto j = i->second.begin(); j != i->second.end(); j++) {
							cout << i->first << ", " << j->first << ", " << j->second.vel << endl;
						}
					}/**/

					/*if (layer != 0)
					{
						for (int x = 0; x < 2; x++)
						{
							for (int y = 0; y < 2; y++)
							{
								vd2d bPos = vd2d(double((fieldx->first << 1 | x) + 0.5), double((fieldy->first << 1 | y) + 0.5)) * (unitDistance >> 1);
								FillCircle((bPos - pos) * zoom + halfScreen, zoom * (unitDistance >> 1) / 2, Pixel((numGravityFields - layer + 1) * (250 / numGravityFields), 0, 0));
							}
						}
					}*/

					if (layer != 0)
					{
						cout << fieldx->first << ", " << fieldy->first << ", " << fieldy->second.vel << " pass to: " << endl;

						for (int x = 0; x < 2; x++)
						{
							findx = layeredGravityFields[layer - 1].find(fieldx->first << 1 + x);

							if (findx != layeredGravityFields[layer - 1].end())
							{
								for (int y = 0; y < 2; y++)
								{
									findy = findx->second.find(fieldy->first << 1 + y);

									if (findy != findx->second.end())
									{
										findy->second.vel += fieldy->second.vel;
										//cout << findx->first << ", " << findy->first << ", " << findy->second.vel << endl;
										/*if (findx->first == 10 && findy->first == 14)
										{
											cout << findy->second.vel << endl;
										}*/
									}
								}
							}
						}
					}
				}
			}
		}

		/*for (auto i = layeredGravityFields[1].begin(); i != layeredGravityFields[1].end(); i++)
		{
			for (auto j = i->second.begin(); j != i->second.end(); j++) {
				cout << i->first << ", " << j->first << ", " << j->second.vel << endl;
			}
		}*/

		for (int i = 0; i < balls.size(); i++)
		{
			cout << balls[i]->pos.touint() << " and " << layeredGravityFields[0][uint32_t(balls[i]->pos.x)][uint32_t(balls[i]->pos.y)].vel << endl;
			balls[i]->vel += layeredGravityFields[0][uint32_t(balls[i]->pos.x)][uint32_t(balls[i]->pos.y)].vel * fElapsedTime;
		}

		for (int i = 0; i < numGravityFields; i++)
		{
			layeredGravityFields[i].clear();
		}
	}

	void ballToBall(ball* ball1, ball* ball2, int dx, int dy, uint32_t partion1x, uint32_t partion1y, uint32_t partion2x, uint32_t partion2y)
	{
		vd2d dpos = vd2d(double(dx) - partion2x + partion1x, double(dy) - partion2y + partion1y) + ball2->pos - ball1->pos;
		double dis = dpos.mag();

		if (dis < ball2->radius + ball1->radius)
		{
			dpos /= dis;
			dis = (ball2->vel - ball1->vel).dot(dpos);

			if (dis < 0)
			{
				dpos *= 2 * dis / (ball1->mass + ball2->mass);
				ball1->vel += dpos * ball2->mass;
				ball2->vel -= dpos * ball1->mass;
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
						ballToBall(partitiony->second[b1], partitiony->second[b2], 0, 0, partitionx->first, partitiony->first, partitionx->first, partitiony->first);

					findy = partitionx->second.find(partitiony->first + 1 & torusMod);

					if (findy != partitionx->second.end())
						for (int b2 = 0; b2 < findy->second.size(); b2++)
							ballToBall(partitiony->second[b1], findy->second[b2], 0, 1, partitionx->first, partitiony->first, partitionx->first, findy->first);

					findx = collisionSpace.find(partitionx->first + 1 & torusMod);

					if (findx != collisionSpace.end())
					{
						findy = findx->second.find(partitiony->first + 1 & torusMod);

						if (findy != findx->second.end())
							for (int b2 = 0; b2 < findy->second.size(); b2++)
								ballToBall(partitiony->second[b1], findy->second[b2], 1, 1, partitionx->first, partitiony->first, findx->first, findy->first);

						findy = findx->second.find(partitiony->first);

						if (findy != findx->second.end())
							for (int b2 = 0; b2 < findy->second.size(); b2++)
								ballToBall(partitiony->second[b1], findy->second[b2], 1, 0, partitionx->first, partitiony->first, findx->first, findy->first);

						findy = findx->second.find(partitiony->first - 1 & torusMod);

						if (findy != findx->second.end())
							for (int b2 = 0; b2 < findy->second.size(); b2++)
								ballToBall(partitiony->second[b1], findy->second[b2], 1, -1, partitionx->first, partitiony->first, findx->first, findy->first);
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
			balls[i]->pos += balls[i]->vel * fElapsedTime;
			balls[i]->pos = balls[i]->pos - (balls[i]->pos / torusRange).floor() * torusRange;
		}
	}

	void drawBalls()
	{
		//Clear(Pixel(0, 0, 0));

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
							vd2d bPos = vd2d(double(partitionx), double(partitiony)) + findy->second[b]->pos - findy->second[b]->pos.touint();

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

		/*for (int i = 0; i < numBalls; i++)
		{
			double randNum = doubleRand() * 6.28318530718;
			vd2d bPos = vd2d(doubleRand(), doubleRand()) * torusRange;
			vd2d bPosv = vd2d(cos(randNum), sin(randNum)) * 1;
			Pixel bColor = mapToRainbow(doubleRand() * 6);
			balls.push_back(new ball(bPos, bPosv, bColor, doubleRand() * (maxRadius - minRadius) + minRadius));
		}*/

		return true;
	}

	bool OnUserUpdate(double fElapsedTime) override
	{
		Clear(Pixel(0, 0, 0));

		fElapsedTime /= 10;

		userControl(fElapsedTime);
		gravity(fElapsedTime);
		collision();
		moveBalls(fElapsedTime);
		drawBalls();

		return true;
	}
};

int main()
{
	/*uint32_t a = 9;

	cout << a - 10 << endl;

	return 0;*/
	Example demo;

	if (demo.Construct(halfScreenx * 2, halfScreeny * 2, 1, 1))
		demo.Start();

	return 0;
}