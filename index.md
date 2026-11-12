---
home: true
title: NGS Course 2026
---

## The case

Three of eleven Weimaraner puppies cannot walk normally: they hop with both hind legs together, like rabbits. Their parents are healthy, and the pathologists found nothing abnormal. During this course you analyse their genomes and find out why, following the study by [Schwarz et al. 2025](https://doi.org/10.1002/mds.30243).

## Exercises

<div class="cards">
{% for d in site.data.days %}
  <a class="card" href="{{ d.url | relative_url }}">
    <span class="card-day">Day {{ d.day }}</span>
    <span class="card-title">{{ d.title }}</span>
    <span class="card-summary">{{ d.summary }}</span>
  </a>
{% endfor %}
</div>

More exercises will be added during the course.
