---
layout: post
title: "The First Stop in Monte Carlo (Methods): Rejection Sampling"
published: true
comments: true
tags: [R, "Monte Carlo"]
---



# Welcome to Monte Carlo

Monte Carlo: site of the Formula 1 Monaco Grand Prix, home to a world-famous casino, and the one of the coolest code names ([coined by John von Neumann](http://en.wikipedia.org/wiki/Monte_Carlo_method#cite_ref-Measure_Anything_pg._46_5-0)) for something mathematical.

According to Wikipedia,

> Monte Carlo methods are a broad class of computational algorithms that rely on repeated random sampling to obtain numerical results.

There are a lot of different applications in fields like finance, engineering, computer graphics, and artificial intelligence. In fact, the best AI for the ancient board game Go relies on Monte Carlo tree searches instead of evaluating all possible moves like most Chess AIs do. You can do a variety of interesting things with Monte Carlo methods, like estimate $$\pi$$ by counting the number samples that fall inside a circle, or make an AI for your game.

In this article, I'm only going to go over one algorithm: **rejection sampling**.

How many times have you been driving down the road thinking, "I really wish I could sample some random variables from a non-standard distribution!" Probably never. That doesn't mean it's not useful. [Here's](http://sumsar.net/blog/2014/10/tiny-data-and-the-socks-of-karl-broman/) a cool example of someone using approximate Bayesian computation (the sampling portion is a form of rejection sampling) to estimate the number of pairs of socks and single socks in a laundry basket based on a very small sample.

So, let's begin!

## The Basics

Let's say we're interested in sampling from some completely off the wall distribution, like say this zig-zag function:

$$
z(x) = \left\{
    \begin{array}{lr}
        x & 0 < x \le 1 \\
        x - 1 & 1 < x \le 2 \\
        0 & \text{otherwise} \\
    \end{array}
    \right.
$$


{% highlight r %}
zig.zag <- function(x) {
    if (0 <= x && x <= 1) x
    else if (1 < x && x <= 2) x - 1
    else 0
}

x <- seq(0, 2, by = 0.01)
y <- sapply(x, zig.zag)

qplot(x, y, geom = 'line')
{% endhighlight %}

![center](/../figs/reject-sample/zig-zag.png)

Does this look even look like a probability distribution function? No, it doesn't. It's not even a continuous function, and with rejection sampling it doesn't even have to be a valid PDF! So how do we go about drawing from it, and how does it work?

### The 10,000ft high perspective

To start off, we draw a box that completely surrounds our zig-zag function. Then, we pick a point at random from our box. If it falls under our zig-zag function, we say we accept this point and return it in the sample. If it falls above the zig-zag function, we reject it and repeat the drawing process. Do this again and again until you get the number of samples that you want. It's a very simple algorithm.

Let's look at an example outcome of these drawings.


{% highlight r %}
# We want 5000 samples from our zig-zag distribution
N <- 5000
i <- 1
draws <- c()
ret <- c()

while (i <= N) {
    # Draw a number along the width of the box
    # (the support of zig-zag)
    T <- runif(1, min=0, max=2)
    # Draw a number along the height of the box
    U <- runif(1, min=0, max=1)

    draws <- rbind(draws, c(T, U))

    # If the point is under the zig-zag,
    # accept it
    if (U <= zig.zag(T)) {
        ret <- c(ret, T)
        i <- i + 1
    }
    # Otherwise, try again
}

summary(ret)
{% endhighlight %}



{% highlight text %}
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
##  0.0261  0.7110  0.9970  1.1700  1.7000  2.0000
{% endhighlight %}

A plot of our data drawn samples and wether or not they were accepted:

![center](/../figs/reject-sample/plot.png)

A histogram of our accepted data points:

![center](/../figs/reject-sample/unnamed-chunk-2.png)

That's pretty close to our distribution. It looks like we're sampling from our zig-zag!

## Down the Rabbit Hole (Diving Into the Details)

Now, I'm going to dive down into the math behind this. I encourage you to keep reading, but if you're just looking for an intuitive feel for this, you can stop here. The above is sufficient for the general idea.

The entire idea of rejection sampling hinges on the idea that if we want to sample from $$f(t)$$ can find an easy to sample density, $$g(t)$$, and a number, $$M$$, such that $$f(t) \le M g(t)$$ for every $$t$$, then for any event, $$E$$,

$$
P(X \in E) = P(T \in E | \text{Accept})
$$

### Great. What's that mean, and can you prove it?

This means that the probability of X being in our event is the same as the probability that T is in our event given that we have accepted this draw.

We'll attempt to prove this by showing they are equivalent.

So if you will all take a minute to remember Bayes' Rule for conditional probability and apply it to our problem:

$$
P(T|Accept) = \frac{P(T\text{ and Accept})}{P(\text{Accept})}
$$

We can now work on simplifying $$P(T\text{ and Accept})$$:

$$
\begin{align}
P(T \text{ and Accept}) &= \int_{E} f_{t, a}(t, a) dt \\
&= \int_{E} f_{A|t}(TRUE | t) g(t) dt \\
&= \int_{E} \frac{f(t)}{C g(t)} g(t) dt \\
&= \frac{1}{C} \int_{E} f(t) dt
\end{align}
$$

So, line by line we have:

1. The definition of a joint probability
2. An application of Bayes' rule
3. Remember, our acceptance algorithm is simply to draw a uniform random variable (pick a random number) from 0 to 1, which we'll call $$U$$, and compare $$U M g(t) \le f(t)$$. This step was just recognizing that the probabilty of accepting given a specific value of t is just the probability that U is greater than the quantity $$\frac{f(t)}{C g(t)}$$. This is a Bernoulli random variable.
4. The last step was just some simplifying algebra.

No we move to the bottom: $$P(\text{Accept})$$

If $$g(T)$$ is defined on the region $$A$$,

$$
\begin{align}
P(\text{Accept}) &= \int_A P_{A,T}(a, t) dt \\
&= \int_A P_{A|T}(TRUE | t) g(t) dt \\
&= \int_A \frac{f(t)}{C g(t)} g(t) dt \\
&= \frac{1}{C} \int_A f(t) dt \\
&= \frac{1}{C}
\end{align}
$$

Again going over this proof line by line:

1. Definition of a marginal probability from a joing PDF
2. Application of Bayes' rule
3. Again the probability of accepting given a specific value of t
4. Simplyfing algebra
5. Integrating over the support of a PDF gives a value of 1

Plugging these back into our original formula leads us to:

$$
\begin{align}
P(X \in E) &= P(T \in E | \text{Accept}) \\
&= P(T|Accept) = \frac{P(T\text{ and Accept})}{P(\text{Accept})} \\
&= \frac{\frac{1}{C} \int_{E} f(t) dt}{\frac{1}{C}} \\
&= \int_E f(t) dt \\
&= P(X \in E)
\end{align}
$$

For that step it's just taking our new pieces, doing some algebra, and applying the definition of the density $$f_X(x)$$ in the last step.

Now that we've proved that those two quantities are equal, we can apply this to any density we want.

### Notes on the choice of a bounding density and scaling constant

The easiest density to sample from is the uniform distribution, so that's what we choose for $$g(t)$$ most of the time. The constant, $$C$$, allows us to scale the height up so that $$C g(t)$$ is greater than any point on the function we want to sample from.

The larger $$C$$ needs to be the larger the proportion of rejected samples will be and the longer the algorithm will take. Increasing $$C$$ increases the amound of rejection area there is. If you have someting that strange like an extremely sharp peak, try bounding with a normal distribution or a gamma distribution with parameters that allow the area between the function to be smaller, e.g. center the normal distribution and give it a small variance to cover the function without having a lot of extra space between.

# Your Turn

Now that you can sample from whatever distribution you want. Can you think of a function you'd want to sample from? What if you wanted to sample from an observed sample distribution? (Those are called empirical density functions.) That's easy now!

As always, if you have any questions or comments, commend below or write me an [e-mail](mailto:thexbarblog@gmail.com).
