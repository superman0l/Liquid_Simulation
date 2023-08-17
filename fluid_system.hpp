#pragma once

#include "fluid_buffer.hpp"
#include "fluid_grid.hpp"

#include <unordered_map>

namespace SPH
{	

	class FluidSystem
	{
	public:
		void init();

		//�������ڹ�ϵ
		void computePressure(void);
		//������ٶ�
		void computeForce(void);
		//�ƶ�
		void advance(void);
		//����Һ��
		void addFluid(Box fluid, float dis);
		//ÿ֡����
		void Update();
		//���������
		void addForce(Eigen::Vector2f a) { m_externalForce += a; };
		void setForce(Eigen::Vector2f a) { m_externalForce = a; };

		//utils
		Buffer getBuffer() { return m_buffer; }
		Box getBound() { return ContainerBox; }

	private:
		GridContainer m_gridContainer;
		Buffer m_buffer;

		//Kernel����
		float m_kernelPoly6;
		float m_kernelSpiky;
		float m_kernelViscosity;

		//ˮ���Ӳ���
		float m_viscosity; //ճ��
		float m_restDensity; //�����ܶ�
		float m_pointMass; //����
		float m_smoothRadius; //�⻬�˰뾶
		float m_gasConstantK; //���峣��K

		float m_speedLimiting; //�ٶ����� 
		float m_accelLimiting; //���ٶ�����
		//�������ڱ߽����
		float m_boundartStiffness;//�����ն�
		float m_boundaryDampening;//�߽�����
		float m_boundary2;//�ӽ������������

		Box FluidBox = Box(Eigen::Vector2f(0.08, 0.04), Eigen::Vector2f(0.12, 0.16)); //Һ��
		Box ContainerBox = Box(Eigen::Vector2f(0.f, 0.f), Eigen::Vector2f(0.2, 0.2)); //�����

		Eigen::Vector2f m_externalForce; //����

	public:
		FluidSystem();
		~FluidSystem();
	};

}