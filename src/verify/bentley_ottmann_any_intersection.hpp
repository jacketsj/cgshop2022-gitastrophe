#pragma once
#include "geom.hpp"
#include <algorithm>
#include <cassert>
#include <exception>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <queue>
#include <random>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace boi {
/**
 * Prepare a vector of segments (in-place),
 * making sure that each segments start point s
 * is to the left (or below) its end point t.
 */
inline void prepare_segments(std::vector<Segment>& segments) {
	for (Segment& s : segments) {
		std::tie(s.s, s.t) = s.left_and_right();
	}
}

/**
 * Draw a random bit with uniform probability.
 * Only uses one RNG call per 64 drawn bits.
 */
inline bool draw_random_bit() {
	struct RngContext {
		RngContext() : rng(std::random_device{}()), bits(0), nbits(0) {}

		bool take() noexcept {
			if (nbits == 0) {
				bits = rng();
				nbits = 64;
			}
			bool res = (bits & 1) != 0;
			bits >>= 1;
			--nbits;
			return res;
		}

		std::mt19937_64 rng;
		std::uint64_t bits;
		unsigned nbits;
	};
	static thread_local RngContext rng;
	return rng.take();
}

/**
 * Implements the Bentley-Ottmann sweep-line algorithm
 * to determine whether a set of segments contains a pair of
 * crossing segments (without approximation).
 * Segments are only counted as intersecting if they share
 * a point that is not an endpoint for at least one of them.
 */
class BentleyOttmannAnyIntersection {
public:
	/**
	 * Information about an intersection (segment indices
	 * and approximate location of the intersection).
	 */
	struct Intersection {
		std::size_t segment_index[2];
		double approximate_location[2];
	};

	/**
	 * Create an instance of the algorithm with a vector of segments.
	 * Does not make a copy of the vector.
	 */
	BentleyOttmannAnyIntersection(const std::vector<Segment>& segments)
			: m_segments(&segments), m_segment_in_list(segments.size(), nullptr),
				m_active_head(nullptr), m_active_tail(nullptr),
				m_intersection(std::nullopt) {}

	~BentleyOttmannAnyIntersection() { skiplist_free(); }

	/**
	 * Perform the actual search for an intersection.
	 */
	std::optional<Intersection> find_intersection() {
		setup_events();
		if (m_active_head) {
			skiplist_free();
		}
		skiplist_init();
		m_intersection = std::nullopt;
		while (!m_events.empty()) {
			Event event = m_events.back();
			m_events.pop_back();
			if (event.entering) {
				handle_entering(event.segment_index);
			} else {
				handle_exiting(event.segment_index);
			}
			if (m_intersection) {
				break;
			}
		}
		return m_intersection;
	}

	/**
	 * Get all intersections (new code)
	 */
	std::vector<Intersection> find_all_intersections() {
		std::vector<Intersection> ret;
		setup_events();
		if (m_active_head) {
			skiplist_free();
		}
		skiplist_init();
		m_intersection = std::nullopt;
		while (!m_events.empty()) {
			Event event = m_events.back();
			m_events.pop_back();
			if (event.entering) {
				handle_entering(event.segment_index);
			} else {
				handle_exiting(event.segment_index);
			}
			if (m_intersection) {
				ret.push_back(m_intersection.value());
			}
		}
		return ret;
	}

private:
	struct SkipNode;

	/**
	 * Check whether the segment given by segment_index intersects the segment
	 * given by node.
	 */
	void check_segment_against(std::size_t segment_index, SkipNode* node) {
		std::size_t other_index = node->segment_index;
		Segment s = (*m_segments)[segment_index], o = (*m_segments)[other_index];
		auto res = do_intersect(s, o);
		if (res) {
			m_intersection =
					Intersection{{segment_index, other_index}, {(*res)[0], (*res)[1]}};
		}
	}

	/**
	 * Handle an event in which a segment enters the set of active segments.
	 */
	void handle_entering(std::size_t segment_index) {
		SkipNode* new_node = skiplist_insert(segment_index);
		if (!new_node)
			return;
		SkipNode *prev = new_node->prev, *next = new_node->next;
		if (prev->prev) {
			check_segment_against(segment_index, prev);
		}
		if (!m_intersection && next->next) {
			check_segment_against(segment_index, next);
		}
	}

	/**
	 * Handle an event in which a segment leaves the set of active segments.
	 */
	void handle_exiting(std::size_t segment_index) {
		SkipNode* rem_node = m_segment_in_list[segment_index];
		SkipNode *prev = rem_node->prev, *next = rem_node->next;
		skiplist_erase(rem_node);
		if (prev->prev && next->next) {
			check_segment_against(prev->segment_index, next);
		}
	}

	/**
	 * Initialize the sweep line data structure that keeps track of the order of
	 * the active segments. In this implementation, we use a skip list for this
	 * purpose.
	 */
	void skiplist_init() {
		std::unique_ptr<SkipNode[]> head{
				SkipNode::allocate(std::numeric_limits<std::size_t>::max(), 1)},
				tail{SkipNode::allocate(std::numeric_limits<std::size_t>::max(), 1)};
		m_active_head = head.release();
		m_active_tail = tail.release();
		m_active_head->next = m_active_tail;
		m_active_tail->prev = m_active_head;
	}

	/**
	 * Free all memory allocated to the sweep line data structure.
	 */
	void skiplist_free() {
		SkipNode* ptr = m_active_head;
		while (ptr) {
			SkipNode* next = ptr->next;
			delete[] ptr;
			ptr = next;
		}
		m_active_head = m_active_tail = nullptr;
	}

	/**
	 * Toss a uniform coin to determine the height of a
	 * newly inserted node in our skip list (between 1 and current_max_level+1).
	 */
	static unsigned skiplist_toss_coins(unsigned current_max_level) {
		unsigned i;
		for (i = 0; i < current_max_level; ++i) {
			if (!draw_random_bit()) {
				break;
			}
		}
		return i + 1;
	}

	/**
	 * Reallocate the skip list head and tail in case the maximum height
	 * increases.
	 */
	void reallocate_head_and_tail(std::size_t new_level) {
		std::unique_ptr<SkipNode[]> new_head{
				SkipNode::allocate(m_active_head->segment_index, new_level)},
				new_tail{SkipNode::allocate(m_active_tail->segment_index, new_level)};
		new_head[new_level - 1].next = &new_tail[new_level - 1];
		new_tail[new_level - 1].prev = &new_head[new_level - 1];
		for (unsigned l = 0; l < new_level - 1; ++l) {
			SkipNode* head = &m_active_head[l];
			SkipNode* tail = &m_active_tail[l];
			SkipNode* nhead = &new_head[l];
			SkipNode* ntail = &new_tail[l];
			SkipNode* next = head->next;
			next->prev = nhead;
			nhead->next = next;
			SkipNode* prev = tail->prev;
			prev->next = ntail;
			ntail->prev = prev;
		}
		delete[] m_active_head;
		delete[] m_active_tail;
		m_active_head = new_head.release();
		m_active_tail = new_tail.release();
	}

	/**
	 * Create a new node for insertion into the skiplist; determines the nodes
	 * height and reallocates head and tail if necessary. The node is not linked
	 * into the skip list.
	 */
	std::unique_ptr<SkipNode[]> create_new_node(std::size_t segment_index) {
		unsigned current_max = m_active_head->num_levels;
		unsigned new_level = skiplist_toss_coins(current_max);
		std::unique_ptr<SkipNode[]> node{
				SkipNode::allocate(segment_index, new_level)};
		if (new_level > current_max) {
			reallocate_head_and_tail(new_level);
		}
		return node;
	}

	/**
	 * Logically negate an optional<bool>, returning an empty optional if the
	 * input is empty.
	 */
	static std::optional<bool> invert(std::optional<bool> input) noexcept {
		if (!input) {
			return input;
		}
		return !*input;
	}

	/**
	 * Set the m_intersection member that stores the intersection found
	 * to the intersection found between segments si1 and si2.
	 */
	void set_intersection(std::size_t si1, std::size_t si2) {
		auto loc =
				approximate_intersection_point((*m_segments)[si1], (*m_segments)[si2]);
		m_intersection = Intersection{{si1, si2}, {loc[0], loc[1]}};
	}

	/**
	 * Compare two segments s1 and s2 for insertion into the sweep line data
	 * structure if segment s1 is vertical. Either of the two segments can be the
	 * newly inserted segment.
	 */
	std::optional<bool> compare_less_at_x_1vert(std::size_t si1, std::size_t si2,
																							Segment s1, Segment s2) {
		if (s2.s.x == s1.s.x) {
			if (s2.s.y <= s1.s.y) {
				return false;
			}
			set_intersection(si1, si2);
			return {};
		}
		std::int64_t delta_x_2 = s2.t.x - s2.s.x;
		std::int64_t delta_y_2 = s2.t.y - s2.s.y;
		std::int64_t delta_x_2_to_now = s1.s.x - s2.s.x;
		std::int64_t lhs = delta_x_2_to_now * delta_y_2;
		std::int64_t rhs_bot = delta_x_2 * (s1.s.y - s2.s.y);
		std::int64_t rhs_top = delta_x_2 * (s1.t.y - s2.s.y);
		if (lhs > rhs_top) {
			return true;
		}
		if (lhs < rhs_bot) {
			return false;
		}
		if (s2.t == s1.t) {
			return true;
		}
		set_intersection(si1, si2);
		return {};
	}

	/**
	 * Compare two segments s1 and s2 for insertion into the sweep line data
	 * structure. s1 is the newly-inserted segment, and neither segment is
	 * vertical.
	 */
	std::optional<bool> compare_less_at_x_novert(std::size_t si1, std::size_t si2,
																							 Segment s1, Segment s2) {
		if (s2.s.x == s1.s.x) {
			if (s1.s.y > s2.s.y) {
				return false;
			}
			// same start point: order by incline
			std::int64_t dx1 = s1.t.x - s1.s.x, dy1 = s1.t.y - s1.s.y;
			std::int64_t dx2 = s2.t.x - s2.s.x, dy2 = s2.t.y - s2.s.y;
			std::int64_t lhs = dy1 * dx2, rhs = dy2 * dx1;
			if (lhs < rhs) {
				return true;
			}
			if (rhs < lhs) {
				return false;
			}
			// overlap
			set_intersection(si1, si2);
			return {};
		}
		std::int64_t delta_x_1_to_now = s1.s.x - s2.s.x;
		std::int64_t delta_x2 = s2.t.x - s2.s.x, delta_y2 = s2.t.y - s2.s.y;
		std::int64_t lhs = delta_x_1_to_now * delta_y2;
		std::int64_t rhs = delta_x2 * (s1.s.y - s2.s.y);
		if (lhs < rhs) {
			return false;
		}
		if (rhs < lhs) {
			return true;
		}
		set_intersection(si1, si2);
		return {};
	}

	/**
	 * Compare two segments for insertion into the sweep line data structure.
	 * si1 is the index of the newly-inserted segment. Returns an empty
	 * optional if the segments cannot be ordered; in that case, an intersection
	 * is always placed into m_intersections.
	 */
	std::optional<bool> compare_less_at_x(std::size_t si1, std::size_t si2) {
		Segment s1 = (*m_segments)[si1], s2 = (*m_segments)[si2];
		bool vert1 = (s1.s.x == s1.t.x), vert2 = (s2.s.x == s2.t.x);
		if (vert1 && vert2) {
			/* since entry events are ordered after exit events, this always is an
			 * overlap */
			set_intersection(si1, si2);
			return {};
		}
		if (vert1) {
			return compare_less_at_x_1vert(si1, si2, s1, s2);
		}
		if (vert2) {
			return invert(compare_less_at_x_1vert(si2, si1, s2, s1));
		}
		return compare_less_at_x_novert(si1, si2, s1, s2);
	}

	/**
	 * Search on the current level of the skip list for the first position where
	 * our newly-inserted segment compares less than the next segment (or the next
	 * segment is the tail). Returns true if an intersection was encountered
	 * during the search. Advances the pointers current/n/nn to the corresponding
	 * position.
	 */
	bool skiplist_search_level(std::size_t segment_index, SkipNode*& current,
														 SkipNode*& n, SkipNode*& nn) {
		for (;;) {
			nn = n->next;
			if (nn == nullptr) {
				return false;
			} else {
				auto res = compare_less_at_x(segment_index, n->segment_index);
				if (!res) {
					return true;
				}
				if (*res) {
					return false;
				}
			}
			current = n;
			n = nn;
		}
	}

	/**
	 * Insert a new segment into the skip list sweep line intersection data
	 * structure. Returns nullptr if an intersection is found during the insertion
	 * process; that intersection is recorded in m_intersection.
	 */
	SkipNode* skiplist_insert(std::size_t segment_index) {
		std::unique_ptr<SkipNode[]> node = create_new_node(segment_index);
		SkipNode* current = &m_active_head[m_active_head->num_levels - 1];
		unsigned l = m_active_head->num_levels - 1;
		unsigned nl = node[0].num_levels - 1;
		SkipNode *n = current->next, *nn;
		// searching above the new nodes' level
		while (l > nl) {
			if (skiplist_search_level(segment_index, current, n, nn)) {
				return nullptr;
			}
			current -= 1;
			n = current->next;
			--l;
		}
		// searching at or below the new nodes' level
		for (;;) {
			if (skiplist_search_level(segment_index, current, n, nn)) {
				return nullptr;
			}
			node[l].prev = current;
			node[l].next = n;
			if (l-- == 0) {
				break;
			}
			current -= 1;
			n = current->next;
		}
		for (l = 0; l < node[0].num_levels; ++l) {
			node[l].prev->next = &node[l];
			node[l].next->prev = &node[l];
		}
		m_segment_in_list[segment_index] = &node[0];
		return node.release();
	}

	/**
	 * Erase a segment from our skip list sweep line data structure.
	 */
	void skiplist_erase(SkipNode* node) {
		for (unsigned i = 0, l = node->num_levels; i < l; ++i) {
			SkipNode* layer = &node[i];
			SkipNode *prev = layer->prev, *next = layer->next;
			prev->next = next;
			next->prev = prev;
		}
		delete[] node;
	}

	/**
	 * Setup the event priority queue.
	 * Implemented as sorted vector of events
	 * since we need not dynamically add events.
	 */
	void setup_events() {
		m_events.clear();
		std::size_t ind = 0;
		for (const auto& seg : *m_segments) {
			m_events.push_back(Event{seg.s, ind, true});
			m_events.push_back(Event{seg.t, ind, false});
			++ind;
		}
		std::sort(m_events.begin(), m_events.end());
	}

	/**
	 * A single node of our skip list sweep line data structure.
	 * A skip list node of height h is allocated as array of h
	 * SkipNodes.
	 */
	struct SkipNode {
		/**
		 * Allocate a new SkipNode of height levels.
		 */
		static SkipNode* allocate(std::size_t segment_index, unsigned levels) {
			SkipNode* nodes = new SkipNode[levels];
			for (unsigned i = 0; i < levels; ++i) {
				nodes[i].segment_index = segment_index;
				nodes[i].level = i;
				nodes[i].num_levels = levels;
				nodes[i].prev = nullptr;
				nodes[i].next = nullptr;
			}
			return nodes;
		}

		std::size_t segment_index;
		SkipNode *prev, *next;
		unsigned level, num_levels;
	};

	/**
	 * Represents an event point, i.e., a segment arriving or leaving.
	 */
	struct Event {
		Event(Point p, std::size_t s, bool e) noexcept
				: point(p), segment_index(s), entering(e) {}

		bool operator<(const Event& other) const noexcept {
			std::int64_t ox = other.point.x, oy = other.point.y, tx = point.x,
									 ty = point.y;
			if (tx < ox) {
				return false;
			} else if (tx > ox) {
				return true;
			}
			if (ty < oy) {
				return false;
			} else if (ty > oy) {
				return true;
			}
			return entering && !other.entering;
		}

		Point point;
		std::size_t segment_index;
		bool entering;
	};

	const std::vector<Segment>* m_segments;
	std::vector<SkipNode*> m_segment_in_list;
	SkipNode *m_active_head, *m_active_tail;
	std::optional<Intersection> m_intersection;
	std::vector<Event> m_events;
};
} // namespace boi
