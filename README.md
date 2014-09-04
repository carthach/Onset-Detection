C implementation of onset detection method using sprectral flux (difference in the spectrum between frames).

Followed a C# tutorial here:-
http://www.badlogicgames.com/wordpress/?p=122

And this paper is also very useful:-
[1] Bello, J. and Daudet, L. 2005. A tutorial on onset detection in music signals. IEEE Transactions on Audio, Speech, and Language Processing. 13, 5 (2005), 1035–1047.

Possibly usable in real-time scenarios (plugins etc.) if you change the peak-picking stage to look back (or delay the incoming signal by a frame or two).
