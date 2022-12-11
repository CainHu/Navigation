#ifndef NAVIGATION_ESKF_UTILS_H
#define NAVIGATION_ESKF_UTILS_H

#include <cmath>
#include <array>
#include <Eigen/Dense>

namespace eskf {
    using namespace std;
    using namespace Eigen;

    template<typename T, unsigned char N>
    class Queue {
    public:
        Queue() = default;;
        explicit Queue(const T &value) {
            for (const T &v : _data) {
                v = value;
            }
        }

        T &operator[](const unsigned char index) {
            return _data[index];
        }

        bool is_empty() const { return _size == 0; };
        bool is_full() const { return _size == N; };
        const T &newest() const { return _data[_head]; };
        const T &oldest() const { return _data[_tail]; };
        unsigned char newest_index() const { return _head; };
        unsigned char oldest_index() const { return _tail; };
        unsigned char size() const { return _size; };
        unsigned char capacity() const { return N; };
        void clear() { _size = 0; _head = N - 1; _tail = 0; }

        void push(const T &value) {
            if (++_head == N) {
                _head = 0;
            }
            _data[_head] = value;

            if (++_size > N) {
                _size = N;
                if (++_tail == N) {
                    _tail = 0;
                }
            }
        }

        bool pop_first_older_than(const unsigned long &timestamp, T & value) {
            for (unsigned char i = 0; i < _size; ++i) {
                unsigned char index = (i > _head) ? (N - i - _head) : _head - i;

                // 离timestamp最接近且延迟小于100ms
                if (timestamp >= _data[index].time_us && timestamp < _data[index].time_us + (unsigned long)1e5) {
                    value = _data[index];

                    // 清空index及index后面的元素
                    _tail = index + 1;
                    if (_tail == N) {
                        _tail = 0;
                    }
                    _size = i;

                    _data[index].time_us = 0;

                    return true;
                }
            }

            // 队列为空 或者 没有timestamp之前的数据, 则返回失败
            return false;
        }

    protected:
        array<T, N> _data {};
        unsigned char _size {0};
        unsigned char _head {N - 1};
        unsigned char _tail {0};
    };

    void rotation_from_axis_angle(Matrix3f &r, const array<float, 3> &a);
    void quaternion_from_axis_angle(Quaternionf &q, const array<float, 3> &a);
    void rotation_from_axis_angle(Matrix3f &r, const Vector3f &a);
    void quaternion_from_axis_angle(Quaternionf &q, const Vector3f &a);
    void normalize_angle(float &angle);
    float normalize_angle(float angle);

}

#endif //NAVIGATION_ESKF_UTILS_H