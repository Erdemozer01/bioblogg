from django import template
from blog.models import Category

register = template.Library()


@register.simple_tag()
def categories_first():
    return Category.objects.all()[:9]


@register.simple_tag()
def categories_middle():
    return Category.objects.all()[9:18]


@register.simple_tag()
def categories_large():
    return Category.objects.all()[18:27]
