#include <array>

namespace eskf {
    using namespace std;

    template<typename T, unsigned char SIZE>
    class Queue {
    public:
        Queue() { };
        Queue(const T& value) {
            for (T & v : _data) {
                v = value;
            }
        };

        T &operator[](const unsigned char index) {
            return _data[index];
        }

        const unsigned char &head() const { return _head; };
        const unsigned char &tail() const { return _tail; };
        const unsigned char &length() const { return SIZE; };

        void push(const T & value) {
            unsigned char head = _head;

            if (!_first_write) {
                head = (_head + 1) % SIZE;
            }

            _data[head] = value;
            _head = head;

            if (_head == _tail && !_first_write) {
                _tail = (_tail + 1) % SIZE;
            } else {
                _first_write = false;
            }
        }

        bool pop_first_older_than(const unsigned long &timestamp, T & value) {
            for (unsigned char i = 0; i < SIZE; ++i) {
                int index = _head - i;
                index = index < 0 ? SIZE + index : index;

                if (timestamp >= _data[index].time_us && timestamp < _data[index].time_us + (unsigned long)1e5) {
                    value = _data[index];

                    // 清空index后面的元素
                    if (index == _head) {   
                        _tail = _head;
                        _first_write = true;
                    } else {    // 清空index
                        _tail = (index + 1) % SIZE;
                    }

                    _data[index].time_us = 0;

                    return true;
                }

                if (index == _tail) {
                    return false;
                }
            }

            return false;
        }

        unsigned char get_length() const {
            return SIZE;
        }

        const T &get_newest() const {
            return _data[_head];
        }

        const T &get_oldest() const {
            return _data[_tail];
        }

        unsigned char get_oldest_index() const {
            return _tail;
        }

        unsigned char get_newest_index() const {
            return _head;
        }

        int entries() const {
            int count = 0;

            for (unsigned char i = 0; i < SIZE; ++i) {
                if (_data[i].time_us != 0) {
                    ++count;
                }
            }

            return count;
        }

    protected:
        bool _first_write {true};
        array<T, SIZE> _data {};
        unsigned char _head {0};
        unsigned char _tail {0};
    };
}